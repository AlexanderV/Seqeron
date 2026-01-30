using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for the GFF/GTF format parser.
/// </summary>
[TestFixture]
public class GffParserTests
{
    #region Test Data

    private const string SimpleGff3 = @"##gff-version 3
chr1	ENSEMBL	gene	1000	5000	.	+	.	ID=gene1;Name=TestGene
chr1	ENSEMBL	mRNA	1000	5000	.	+	.	ID=transcript1;Parent=gene1
chr1	ENSEMBL	exon	1000	1500	.	+	.	ID=exon1;Parent=transcript1
chr1	ENSEMBL	exon	2000	2500	.	+	.	ID=exon2;Parent=transcript1
chr1	ENSEMBL	CDS	1100	1500	.	+	0	ID=cds1;Parent=transcript1
chr1	ENSEMBL	CDS	2000	2400	.	+	2	ID=cds2;Parent=transcript1";

    private const string SimpleGtf = "chr1\tENSEMBL\tgene\t1000\t5000\t.\t+\t.\tgene_id \"ENSG00001\"; gene_name \"TestGene\";\n" +
        "chr1\tENSEMBL\ttranscript\t1000\t5000\t.\t+\t.\tgene_id \"ENSG00001\"; transcript_id \"ENST00001\";\n" +
        "chr1\tENSEMBL\texon\t1000\t1500\t.\t+\t.\tgene_id \"ENSG00001\"; transcript_id \"ENST00001\"; exon_number \"1\";\n" +
        "chr1\tENSEMBL\texon\t2000\t2500\t.\t+\t.\tgene_id \"ENSG00001\"; transcript_id \"ENST00001\"; exon_number \"2\";";

    private const string MultiChromGff = @"##gff-version 3
chr1	.	gene	1000	2000	.	+	.	ID=gene1
chr1	.	gene	3000	4000	.	-	.	ID=gene2
chr2	.	gene	1000	2000	.	+	.	ID=gene3
chr2	.	gene	5000	6000	.	+	.	ID=gene4";

    #endregion

    #region Basic Parsing Tests

    [Test]
    public void Parse_GFF3_ReturnsCorrectRecords()
    {
        var records = GffParser.Parse(SimpleGff3, GffParser.GffFormat.GFF3).ToList();

        Assert.That(records, Has.Count.EqualTo(6));
        Assert.That(records[0].Type, Is.EqualTo("gene"));
        Assert.That(records[0].Start, Is.EqualTo(1000));
        Assert.That(records[0].End, Is.EqualTo(5000));
    }

    [Test]
    public void Parse_GTF_ReturnsCorrectRecords()
    {
        var records = GffParser.Parse(SimpleGtf, GffParser.GffFormat.GTF).ToList();

        Assert.That(records, Has.Count.GreaterThanOrEqualTo(4));
    }

    [Test]
    public void Parse_EmptyContent_ReturnsEmpty()
    {
        var records = GffParser.Parse("", GffParser.GffFormat.Auto).ToList();
        Assert.That(records, Is.Empty);
    }

    [Test]
    public void Parse_NullContent_ReturnsEmpty()
    {
        var records = GffParser.Parse((string)null!, GffParser.GffFormat.Auto).ToList();
        Assert.That(records, Is.Empty);
    }

    [Test]
    public void Parse_SkipsComments()
    {
        const string gff = @"##gff-version 3
# This is a comment
chr1	.	gene	1000	2000	.	+	.	ID=gene1";

        var records = GffParser.Parse(gff).ToList();
        Assert.That(records, Has.Count.EqualTo(1));
    }

    [Test]
    public void Parse_SkipsEmptyLines()
    {
        const string gff = @"##gff-version 3
chr1	.	gene	1000	2000	.	+	.	ID=gene1

chr1	.	gene	3000	4000	.	+	.	ID=gene2";

        var records = GffParser.Parse(gff).ToList();
        Assert.That(records, Has.Count.EqualTo(2));
    }

    #endregion

    #region Attribute Parsing Tests

    [Test]
    public void Parse_GFF3Attributes_ParsedCorrectly()
    {
        var records = GffParser.Parse(SimpleGff3, GffParser.GffFormat.GFF3).ToList();
        var gene = records.First(r => r.Type == "gene");

        Assert.That(gene.Attributes.ContainsKey("ID"), Is.True);
        Assert.That(gene.Attributes["ID"], Is.EqualTo("gene1"));
        Assert.That(gene.Attributes["Name"], Is.EqualTo("TestGene"));
    }

    [Test]
    public void Parse_GTFAttributes_ParsedCorrectly()
    {
        var records = GffParser.Parse(SimpleGtf, GffParser.GffFormat.GTF).ToList();
        var gene = records.First(r => r.Type == "gene");

        Assert.That(gene.Attributes.ContainsKey("gene_id"), Is.True);
        Assert.That(gene.Attributes["gene_id"], Is.EqualTo("ENSG00001"));
    }

    [Test]
    public void GetAttribute_ExistingAttribute_ReturnsValue()
    {
        var records = GffParser.Parse(SimpleGff3).ToList();
        var gene = records.First(r => r.Type == "gene");

        var name = GffParser.GetAttribute(gene, "Name");
        Assert.That(name, Is.EqualTo("TestGene"));
    }

    [Test]
    public void GetAttribute_NonexistentAttribute_ReturnsNull()
    {
        var records = GffParser.Parse(SimpleGff3).ToList();
        var gene = records.First(r => r.Type == "gene");

        var value = GffParser.GetAttribute(gene, "NonExistent");
        Assert.That(value, Is.Null);
    }

    [Test]
    public void GetGeneName_ReturnsGeneName()
    {
        var records = GffParser.Parse(SimpleGff3).ToList();
        var gene = records.First(r => r.Type == "gene");

        var name = GffParser.GetGeneName(gene);
        Assert.That(name, Is.EqualTo("TestGene"));
    }

    #endregion

    #region Filtering Tests

    [Test]
    public void FilterByType_ReturnsMatchingTypes()
    {
        var records = GffParser.Parse(SimpleGff3).ToList();
        var genes = GffParser.FilterByType(records, "gene").ToList();

        Assert.That(genes, Has.Count.EqualTo(1));
        Assert.That(genes.All(r => r.Type == "gene"), Is.True);
    }

    [Test]
    public void FilterByType_MultipleTypes_ReturnsAll()
    {
        var records = GffParser.Parse(SimpleGff3).ToList();
        var filtered = GffParser.FilterByType(records, "gene", "exon").ToList();

        Assert.That(filtered, Has.Count.EqualTo(3)); // 1 gene + 2 exons
    }

    [Test]
    public void FilterBySeqid_ReturnsMatchingChromosome()
    {
        var records = GffParser.Parse(MultiChromGff).ToList();
        var chr1 = GffParser.FilterBySeqid(records, "chr1").ToList();

        Assert.That(chr1, Has.Count.EqualTo(2));
        Assert.That(chr1.All(r => r.Seqid == "chr1"), Is.True);
    }

    [Test]
    public void FilterByRegion_ReturnsOverlappingFeatures()
    {
        var records = GffParser.Parse(MultiChromGff).ToList();
        var inRegion = GffParser.FilterByRegion(records, "chr1", 1500, 3500).ToList();

        Assert.That(inRegion, Has.Count.EqualTo(2)); // gene1 and gene2 overlap this region
    }

    [Test]
    public void GetGenes_ReturnsOnlyGenes()
    {
        var records = GffParser.Parse(SimpleGff3).ToList();
        var genes = GffParser.GetGenes(records).ToList();

        Assert.That(genes, Has.Count.EqualTo(1));
        Assert.That(genes[0].Type, Is.EqualTo("gene"));
    }

    [Test]
    public void GetExons_ReturnsOnlyExons()
    {
        var records = GffParser.Parse(SimpleGff3).ToList();
        var exons = GffParser.GetExons(records).ToList();

        Assert.That(exons, Has.Count.EqualTo(2));
        Assert.That(exons.All(e => e.Type == "exon"), Is.True);
    }

    [Test]
    public void GetCDS_ReturnsOnlyCDS()
    {
        var records = GffParser.Parse(SimpleGff3).ToList();
        var cds = GffParser.GetCDS(records).ToList();

        Assert.That(cds, Has.Count.EqualTo(2));
        Assert.That(cds.All(c => c.Type == "CDS"), Is.True);
    }

    #endregion

    #region Gene Model Building Tests

    [Test]
    public void BuildGeneModels_CreatesHierarchy()
    {
        var records = GffParser.Parse(SimpleGff3).ToList();
        var models = GffParser.BuildGeneModels(records).ToList();

        Assert.That(models, Has.Count.EqualTo(1));
        Assert.That(models[0].Gene.Type, Is.EqualTo("gene"));
        Assert.That(models[0].Transcripts.Count, Is.GreaterThanOrEqualTo(1));
        Assert.That(models[0].Exons.Count, Is.GreaterThanOrEqualTo(2));
        Assert.That(models[0].CDS.Count, Is.GreaterThanOrEqualTo(2));
    }

    [Test]
    public void BuildGeneModels_MultipleGenes_BuildsAll()
    {
        var records = GffParser.Parse(MultiChromGff).ToList();
        var models = GffParser.BuildGeneModels(records).ToList();

        Assert.That(models, Has.Count.EqualTo(4));
    }

    #endregion

    #region Statistics Tests

    [Test]
    public void CalculateStatistics_ReturnsCorrectCounts()
    {
        var records = GffParser.Parse(SimpleGff3).ToList();
        var stats = GffParser.CalculateStatistics(records);

        Assert.That(stats.TotalFeatures, Is.EqualTo(6));
        Assert.That(stats.GeneCount, Is.EqualTo(1));
        Assert.That(stats.ExonCount, Is.EqualTo(2));
    }

    [Test]
    public void CalculateStatistics_FeatureTypeCounts_Correct()
    {
        var records = GffParser.Parse(SimpleGff3).ToList();
        var stats = GffParser.CalculateStatistics(records);

        Assert.That(stats.FeatureTypeCounts["gene"], Is.EqualTo(1));
        Assert.That(stats.FeatureTypeCounts["exon"], Is.EqualTo(2));
        Assert.That(stats.FeatureTypeCounts["CDS"], Is.EqualTo(2));
    }

    [Test]
    public void CalculateStatistics_SequenceIds_Listed()
    {
        var records = GffParser.Parse(MultiChromGff).ToList();
        var stats = GffParser.CalculateStatistics(records);

        Assert.That(stats.SequenceIds, Contains.Item("chr1"));
        Assert.That(stats.SequenceIds, Contains.Item("chr2"));
    }

    #endregion

    #region Writing Tests

    [Test]
    public void WriteToStream_GFF3Format_ValidOutput()
    {
        var records = GffParser.Parse(SimpleGff3).ToList();
        using var writer = new StringWriter();

        GffParser.WriteToStream(writer, records, GffParser.GffFormat.GFF3);
        var output = writer.ToString();

        Assert.That(output, Does.Contain("##gff-version 3"));
        Assert.That(output, Does.Contain("chr1"));
        Assert.That(output, Does.Contain("gene"));
    }

    [Test]
    public void WriteAndRead_Roundtrip_PreservesData()
    {
        var original = GffParser.Parse(SimpleGff3).ToList();
        using var writer = new StringWriter();

        GffParser.WriteToStream(writer, original, GffParser.GffFormat.GFF3);
        var output = writer.ToString();

        var parsed = GffParser.Parse(output, GffParser.GffFormat.GFF3).ToList();

        Assert.That(parsed.Count, Is.EqualTo(original.Count));
        Assert.That(parsed[0].Seqid, Is.EqualTo(original[0].Seqid));
        Assert.That(parsed[0].Start, Is.EqualTo(original[0].Start));
        Assert.That(parsed[0].End, Is.EqualTo(original[0].End));
    }

    #endregion

    #region Utility Tests

    [Test]
    public void ExtractSequence_PlusStrand_ReturnsSequence()
    {
        var record = new GffParser.GffRecord(
            "chr1", ".", "exon", 5, 10, null, '+', null, new Dictionary<string, string>());

        var reference = "ACGTACGTACGTACGT";
        var sequence = GffParser.ExtractSequence(record, reference);

        // GFF is 1-based: positions 5-10 → C# indices [4..10] → "ACGTAC" (6 chars)
        Assert.That(sequence, Is.EqualTo("ACGTAC"));
    }

    [Test]
    public void ExtractSequence_MinusStrand_ReturnsReverseComplement()
    {
        var record = new GffParser.GffRecord(
            "chr1", ".", "exon", 1, 4, null, '-', null, new Dictionary<string, string>());

        var reference = "ACGT";
        var sequence = GffParser.ExtractSequence(record, reference);

        Assert.That(sequence, Is.EqualTo("ACGT")); // Reverse complement of ACGT
    }

    [Test]
    public void MergeOverlapping_MergesCorrectly()
    {
        const string gff = @"##gff-version 3
chr1	.	exon	100	200	.	+	.	ID=e1
chr1	.	exon	150	250	.	+	.	ID=e2
chr1	.	exon	400	500	.	+	.	ID=e3";

        var records = GffParser.Parse(gff).ToList();
        var merged = GffParser.MergeOverlapping(records).ToList();

        Assert.That(merged, Has.Count.EqualTo(2));
        Assert.That(merged[0].Start, Is.EqualTo(100));
        Assert.That(merged[0].End, Is.EqualTo(250)); // Merged
        Assert.That(merged[1].Start, Is.EqualTo(400));
    }

    #endregion

    #region Edge Cases

    [Test]
    public void Parse_MalformedLine_Skips()
    {
        const string gff = @"##gff-version 3
chr1	.	gene
chr1	.	gene	1000	2000	.	+	.	ID=gene1";

        var records = GffParser.Parse(gff).ToList();
        Assert.That(records, Has.Count.EqualTo(1));
    }

    [Test]
    public void Parse_SpecialCharacters_Unescaped()
    {
        const string gff = @"##gff-version 3
chr1	.	gene	1000	2000	.	+	.	ID=gene1;Name=Test%3BGene";

        var records = GffParser.Parse(gff).ToList();
        // URL-encoded ; should be unescaped
        Assert.That(records[0].Attributes["Name"], Does.Contain("Gene"));
    }

    [Test]
    public void Parse_NoScore_ScoreIsNull()
    {
        var records = GffParser.Parse(SimpleGff3).ToList();
        Assert.That(records[0].Score, Is.Null);
    }

    [Test]
    public void Parse_WithScore_ScoreIsParsed()
    {
        const string gff = @"chr1	.	gene	1000	2000	99.5	+	.	ID=gene1";

        var records = GffParser.Parse(gff).ToList();
        Assert.That(records[0].Score, Is.EqualTo(99.5).Within(0.01));
    }

    #endregion

    #region File I/O Tests

    [Test]
    public void ParseFile_NonexistentFile_ReturnsEmpty()
    {
        var records = GffParser.ParseFile("nonexistent.gff").ToList();
        Assert.That(records, Is.Empty);
    }

    [Test]
    public void ParseFile_ValidFile_ParsesRecords()
    {
        var tempFile = Path.GetTempFileName();
        try
        {
            File.WriteAllText(tempFile, SimpleGff3);
            var records = GffParser.ParseFile(tempFile).ToList();

            Assert.That(records, Has.Count.EqualTo(6));
        }
        finally
        {
            File.Delete(tempFile);
        }
    }

    #endregion
}
