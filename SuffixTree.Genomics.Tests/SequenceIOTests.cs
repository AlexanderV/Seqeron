using NUnit.Framework;
using SuffixTree.Genomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

[TestFixture]
public class SequenceIOTests
{
    #region GenBank Tests

    private const string SampleGenBank = @"LOCUS       TEST001               100 bp    DNA     linear   UNK 01-JAN-2024
DEFINITION  Test sequence for unit testing.
ACCESSION   TEST001
SOURCE      Test organism
  ORGANISM  Test organism
            Bacteria; Proteobacteria.
FEATURES             Location/Qualifiers
     gene            1..50
                     /gene=""testGene""
                     /locus_tag=""TEST_001""
     CDS             10..40
                     /gene=""testGene""
                     /product=""test protein""
     gene            complement(60..90)
                     /gene=""reverseGene""
ORIGIN
        1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat
       51 gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc
//
";

    [Test]
    public void ParseGenBankString_ValidRecord_ParsesCorrectly()
    {
        var records = SequenceIO.ParseGenBankString(SampleGenBank).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        Assert.That(records[0].Id, Is.EqualTo("TEST001"));
        Assert.That(records[0].Accession, Is.EqualTo("TEST001"));
        Assert.That(records[0].Description, Is.EqualTo("Test sequence for unit testing."));
        Assert.That(records[0].Organism, Is.EqualTo("Test organism"));
    }

    [Test]
    public void ParseGenBankString_ValidRecord_ParsesSequence()
    {
        var records = SequenceIO.ParseGenBankString(SampleGenBank).ToList();

        Assert.That(records[0].Sequence.Length, Is.EqualTo(100));
        Assert.That(records[0].Sequence, Does.StartWith("ATGCATGCAT"));
    }

    [Test]
    public void ParseGenBankString_ValidRecord_ParsesFeatures()
    {
        var records = SequenceIO.ParseGenBankString(SampleGenBank).ToList();

        Assert.That(records[0].Features, Has.Count.EqualTo(3));

        var gene = records[0].Features.First(f => f.Type == "gene" && f.Strand == '+');
        Assert.That(gene.Start, Is.EqualTo(1));
        Assert.That(gene.End, Is.EqualTo(50));

        var cds = records[0].Features.First(f => f.Type == "CDS");
        Assert.That(cds.Qualifiers["product"], Is.EqualTo("test protein"));
    }

    [Test]
    public void ParseGenBankString_ComplementFeature_ParsesStrand()
    {
        var records = SequenceIO.ParseGenBankString(SampleGenBank).ToList();

        var reverseGene = records[0].Features.First(f => f.Qualifiers.ContainsKey("gene") &&
                                                          f.Qualifiers["gene"] == "reverseGene");
        Assert.That(reverseGene.Strand, Is.EqualTo('-'));
    }

    [Test]
    public void ParseGenBankString_EmptyString_ReturnsEmpty()
    {
        var records = SequenceIO.ParseGenBankString("").ToList();
        Assert.That(records, Is.Empty);
    }

    [Test]
    public void ToGenBank_ValidRecord_GeneratesValidFormat()
    {
        var record = new SequenceIO.SequenceRecord(
            Id: "SEQ001",
            Accession: "ACC001",
            Description: "Test sequence",
            Sequence: "ATGCATGCATGC",
            Organism: "Test organism",
            Taxonomy: "Test; Taxonomy",
            Date: new DateTime(2024, 1, 15),
            Features: new List<SequenceIO.SequenceFeature>
            {
                new("gene", 1, 12, '+', new Dictionary<string, string> { { "gene", "test" } })
            },
            References: new List<SequenceIO.Reference>(),
            Metadata: new Dictionary<string, string>());

        string genbank = SequenceIO.ToGenBank(record);

        Assert.That(genbank, Does.Contain("LOCUS"));
        Assert.That(genbank, Does.Contain("SEQ001"));
        Assert.That(genbank, Does.Contain("ORIGIN"));
        Assert.That(genbank, Does.Contain("atgcatgcat"));
        Assert.That(genbank, Does.Contain("//"));
    }

    [Test]
    public void ParseGenBank_RoundTrip_PreservesData()
    {
        var original = new SequenceIO.SequenceRecord(
            Id: "ROUND001",
            Accession: "ACC001",
            Description: "Round trip test",
            Sequence: "ATGCATGCATGCATGCATGCATGCATGCATGC",
            Organism: "Test",
            Taxonomy: null,
            Date: null,
            Features: new List<SequenceIO.SequenceFeature>(),
            References: new List<SequenceIO.Reference>(),
            Metadata: new Dictionary<string, string>());

        string genbank = SequenceIO.ToGenBank(original);
        var parsed = SequenceIO.ParseGenBankString(genbank).First();

        Assert.That(parsed.Id, Is.EqualTo("ROUND001"));
        Assert.That(parsed.Sequence, Is.EqualTo(original.Sequence));
    }

    #endregion

    #region EMBL Tests

    private const string SampleEmbl = @"ID   TEST001; SV 1; linear; DNA; STD; UNK; 100 BP.
AC   AC123456;
DE   Test sequence description.
OS   Test organism
SQ   Sequence 100 BP;
     atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat        50
     gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc       100
//
";

    [Test]
    public void ParseEmblString_ValidRecord_ParsesCorrectly()
    {
        var records = SequenceIO.ParseEmblString(SampleEmbl).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        Assert.That(records[0].Id, Is.EqualTo("TEST001"));
        Assert.That(records[0].Accession, Is.EqualTo("AC123456"));
        Assert.That(records[0].Organism, Is.EqualTo("Test organism"));
    }

    [Test]
    public void ParseEmblString_ValidRecord_ParsesSequence()
    {
        var records = SequenceIO.ParseEmblString(SampleEmbl).ToList();

        Assert.That(records[0].Sequence.Length, Is.EqualTo(100));
        Assert.That(records[0].Sequence, Does.StartWith("ATGCATGCAT"));
    }

    [Test]
    public void ParseEmblString_EmptyString_ReturnsEmpty()
    {
        var records = SequenceIO.ParseEmblString("").ToList();
        Assert.That(records, Is.Empty);
    }

    #endregion

    #region BED Tests

    private const string SampleBed = @"# Comment line
chr1	100	200	gene1	500	+	110	190	255,0,0
chr1	300	400	gene2	600	-
chr2	500	600
";

    [Test]
    public void ParseBedString_ValidEntries_ParsesCorrectly()
    {
        var entries = SequenceIO.ParseBedString(SampleBed).ToList();

        Assert.That(entries, Has.Count.EqualTo(3));
        Assert.That(entries[0].Chromosome, Is.EqualTo("chr1"));
        Assert.That(entries[0].Start, Is.EqualTo(100));
        Assert.That(entries[0].End, Is.EqualTo(200));
        Assert.That(entries[0].Name, Is.EqualTo("gene1"));
        Assert.That(entries[0].Score, Is.EqualTo(500));
        Assert.That(entries[0].Strand, Is.EqualTo('+'));
    }

    [Test]
    public void ParseBedString_MinimalEntry_ParsesCorrectly()
    {
        var entries = SequenceIO.ParseBedString(SampleBed).ToList();

        Assert.That(entries[2].Chromosome, Is.EqualTo("chr2"));
        Assert.That(entries[2].Start, Is.EqualTo(500));
        Assert.That(entries[2].End, Is.EqualTo(600));
        Assert.That(entries[2].Name, Is.Null);
    }

    [Test]
    public void ParseBedString_SkipsComments()
    {
        string bed = "# Comment\nchr1\t100\t200\n";
        var entries = SequenceIO.ParseBedString(bed).ToList();

        Assert.That(entries, Has.Count.EqualTo(1));
    }

    [Test]
    public void ToBed_ValidEntries_GeneratesCorrectFormat()
    {
        var entries = new List<SequenceIO.BedEntry>
        {
            new("chr1", 100, 200, "gene1", 500, '+', 110, 190, "255,0,0")
        };

        string bed = SequenceIO.ToBed(entries);

        Assert.That(bed, Does.Contain("chr1\t100\t200\tgene1"));
    }

    [Test]
    public void ParseBed_RoundTrip_PreservesData()
    {
        var original = new List<SequenceIO.BedEntry>
        {
            new("chr1", 100, 200, "gene1", 500, '+')
        };

        string bed = SequenceIO.ToBed(original);
        var parsed = SequenceIO.ParseBedString(bed).ToList();

        Assert.That(parsed[0].Chromosome, Is.EqualTo("chr1"));
        Assert.That(parsed[0].Start, Is.EqualTo(100));
        Assert.That(parsed[0].End, Is.EqualTo(200));
    }

    #endregion

    #region GFF Tests

    private const string SampleGff = @"##gff-version 3
chr1	source	gene	100	500	.	+	.	ID=gene1;Name=TestGene
chr1	source	CDS	150	450	100.5	+	0	Parent=gene1;product=TestProtein
chr1	source	exon	100	200	.	-	.	ID=exon1
";

    [Test]
    public void ParseGffString_ValidEntries_ParsesCorrectly()
    {
        var entries = SequenceIO.ParseGffString(SampleGff).ToList();

        Assert.That(entries, Has.Count.EqualTo(3));
        Assert.That(entries[0].SeqId, Is.EqualTo("chr1"));
        Assert.That(entries[0].Type, Is.EqualTo("gene"));
        Assert.That(entries[0].Start, Is.EqualTo(100));
        Assert.That(entries[0].End, Is.EqualTo(500));
        Assert.That(entries[0].Strand, Is.EqualTo('+'));
    }

    [Test]
    public void ParseGffString_ParsesAttributes()
    {
        var entries = SequenceIO.ParseGffString(SampleGff).ToList();

        Assert.That(entries[0].Attributes["ID"], Is.EqualTo("gene1"));
        Assert.That(entries[0].Attributes["Name"], Is.EqualTo("TestGene"));
    }

    [Test]
    public void ParseGffString_ParsesScore()
    {
        var entries = SequenceIO.ParseGffString(SampleGff).ToList();

        Assert.That(entries[0].Score, Is.Null);
        Assert.That(entries[1].Score, Is.EqualTo(100.5).Within(0.1));
    }

    [Test]
    public void ParseGffString_ParsesPhase()
    {
        var entries = SequenceIO.ParseGffString(SampleGff).ToList();

        Assert.That(entries[0].Phase, Is.Null);
        Assert.That(entries[1].Phase, Is.EqualTo(0));
    }

    [Test]
    public void ToGff_ValidEntries_GeneratesGff3Format()
    {
        var entries = new List<SequenceIO.GffEntry>
        {
            new("chr1", "src", "gene", 100, 500, null, '+', null,
                new Dictionary<string, string> { { "ID", "gene1" } })
        };

        string gff = SequenceIO.ToGff(entries, gff3: true);

        Assert.That(gff, Does.Contain("##gff-version 3"));
        Assert.That(gff, Does.Contain("chr1\tsrc\tgene\t100\t500"));
    }

    [Test]
    public void ParseGff_SkipsComments()
    {
        string gff = "##gff-version 3\n# Comment\nchr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1\n";
        var entries = SequenceIO.ParseGffString(gff).ToList();

        Assert.That(entries, Has.Count.EqualTo(1));
    }

    #endregion

    #region SAM Tests

    private const string SampleSam = @"@HD	VN:1.0	SO:coordinate
@SQ	SN:chr1	LN:1000000
read1	99	chr1	100	60	50M	=	200	150	ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2	163	chr1	200	40	30M5I15M	=	100	-150	ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
";

    [Test]
    public void ParseSamString_ValidRecords_ParsesCorrectly()
    {
        var records = SequenceIO.ParseSamString(SampleSam).ToList();

        Assert.That(records, Has.Count.EqualTo(2));
        Assert.That(records[0].ReadName, Is.EqualTo("read1"));
        Assert.That(records[0].Position, Is.EqualTo(100));
        Assert.That(records[0].MappingQuality, Is.EqualTo(60));
        Assert.That(records[0].Cigar, Is.EqualTo("50M"));
    }

    [Test]
    public void ParseSamString_SkipsHeaders()
    {
        var records = SequenceIO.ParseSamString(SampleSam).ToList();
        Assert.That(records, Has.Count.EqualTo(2)); // Only alignments, not headers
    }

    [Test]
    public void ParseSamString_ParsesFlags()
    {
        var records = SequenceIO.ParseSamString(SampleSam).ToList();

        Assert.That(records[0].Flag, Is.EqualTo(99));
        Assert.That(SequenceIO.IsPaired(99), Is.True);
        Assert.That(SequenceIO.IsProperPair(99), Is.True);
    }

    [Test]
    public void SamFlagChecks_WorkCorrectly()
    {
        Assert.That(SequenceIO.IsPaired(1), Is.True);
        Assert.That(SequenceIO.IsPaired(0), Is.False);
        Assert.That(SequenceIO.IsUnmapped(4), Is.True);
        Assert.That(SequenceIO.IsReverse(16), Is.True);
        Assert.That(SequenceIO.IsRead1(64), Is.True);
        Assert.That(SequenceIO.IsRead2(128), Is.True);
        Assert.That(SequenceIO.IsSecondary(256), Is.True);
        Assert.That(SequenceIO.IsDuplicate(1024), Is.True);
        Assert.That(SequenceIO.IsSupplementary(2048), Is.True);
    }

    #endregion

    #region VCF Tests

    private const string SampleVcf = @"##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	100	rs123	A	G	30.5	PASS	DP=50;AF=0.5	GT:DP	0/1:25
chr1	200	.	AT	A	.	.	DP=30	GT	1/1
chr2	300	rs456	C	T,G	45.0	PASS	DP=100;MQ=60	GT	0/1
";

    [Test]
    public void ParseVcfString_ValidRecords_ParsesCorrectly()
    {
        var records = SequenceIO.ParseVcfString(SampleVcf).ToList();

        Assert.That(records, Has.Count.EqualTo(3));
        Assert.That(records[0].Chromosome, Is.EqualTo("chr1"));
        Assert.That(records[0].Position, Is.EqualTo(100));
        Assert.That(records[0].Reference, Is.EqualTo("A"));
        Assert.That(records[0].Alternatives[0], Is.EqualTo("G"));
    }

    [Test]
    public void ParseVcfString_ParsesQuality()
    {
        var records = SequenceIO.ParseVcfString(SampleVcf).ToList();

        Assert.That(records[0].Quality, Is.EqualTo(30.5).Within(0.1));
        Assert.That(records[1].Quality, Is.Null);
    }

    [Test]
    public void ParseVcfString_ParsesInfo()
    {
        var records = SequenceIO.ParseVcfString(SampleVcf).ToList();

        Assert.That(records[0].Info["DP"], Is.EqualTo("50"));
        Assert.That(records[0].Info["AF"], Is.EqualTo("0.5"));
    }

    [Test]
    public void ParseVcfString_MultipleAlternatives_ParsesAll()
    {
        var records = SequenceIO.ParseVcfString(SampleVcf).ToList();

        Assert.That(records[2].Alternatives, Has.Count.EqualTo(2));
        Assert.That(records[2].Alternatives[0], Is.EqualTo("T"));
        Assert.That(records[2].Alternatives[1], Is.EqualTo("G"));
    }

    [Test]
    public void ParseVcfString_SkipsHeaders()
    {
        var records = SequenceIO.ParseVcfString(SampleVcf).ToList();
        Assert.That(records, Has.Count.EqualTo(3));
    }

    #endregion

    #region Phylip Tests

    private const string SamplePhylip = @" 3 20
Seq1      ATGCATGCATGCATGCATGC
Seq2      ATGCATGCATGCATGCATGC
Seq3      TTTTTTTTTTTTTTTTTTTT
";

    [Test]
    public void ParsePhylip_ValidFile_ParsesCorrectly()
    {
        using var reader = new StringReader(SamplePhylip);
        var sequences = SequenceIO.ParsePhylip(reader).ToList();

        Assert.That(sequences, Has.Count.EqualTo(3));
        Assert.That(sequences[0].Name, Is.EqualTo("Seq1"));
        Assert.That(sequences[0].Sequence, Is.EqualTo("ATGCATGCATGCATGCATGC"));
    }

    [Test]
    public void ToPhylip_ValidSequences_GeneratesCorrectFormat()
    {
        var sequences = new List<(string Name, string Sequence)>
        {
            ("Seq1", "ATGCATGC"),
            ("Seq2", "GCTAGCTA")
        };

        string phylip = SequenceIO.ToPhylip(sequences);

        Assert.That(phylip, Does.Contain("2 8"));
        Assert.That(phylip, Does.Contain("Seq1"));
    }

    [Test]
    public void ToPhylip_Interleaved_GeneratesCorrectFormat()
    {
        var sequences = new List<(string Name, string Sequence)>
        {
            ("Seq1", string.Concat(Enumerable.Repeat("ATGC", 20))),
            ("Seq2", string.Concat(Enumerable.Repeat("GCTA", 20)))
        };

        string phylip = SequenceIO.ToPhylip(sequences, interleaved: true);

        Assert.That(phylip, Does.Contain("2 80"));
    }

    [Test]
    public void ToPhylip_EmptySequences_ReturnsEmpty()
    {
        var sequences = new List<(string Name, string Sequence)>();
        string phylip = SequenceIO.ToPhylip(sequences);

        Assert.That(phylip, Is.EqualTo(""));
    }

    #endregion

    #region Clustal Tests

    private const string SampleClustal = @"CLUSTAL W (1.83) multiple sequence alignment

Seq1      ATGCATGCAT
Seq2      ATGCATGCAT
Seq3      TTTTTTTTTT
          **** *****

Seq1      ATGCATGCAT
Seq2      ATGCATGCAT
Seq3      TTTTTTTTTT
";

    [Test]
    public void ParseClustal_ValidFile_ParsesCorrectly()
    {
        using var reader = new StringReader(SampleClustal);
        var sequences = SequenceIO.ParseClustal(reader).ToList();

        Assert.That(sequences, Has.Count.EqualTo(3));
        Assert.That(sequences[0].Name, Is.EqualTo("Seq1"));
        Assert.That(sequences[0].Sequence, Is.EqualTo("ATGCATGCATATGCATGCAT")); // Combined blocks
    }

    [Test]
    public void ParseClustal_SkipsHeader()
    {
        using var reader = new StringReader(SampleClustal);
        var sequences = SequenceIO.ParseClustal(reader).ToList();

        Assert.That(sequences.All(s => !s.Name.StartsWith("CLUSTAL")));
    }

    [Test]
    public void ParseClustal_SkipsConservationLines()
    {
        using var reader = new StringReader(SampleClustal);
        var sequences = SequenceIO.ParseClustal(reader).ToList();

        Assert.That(sequences.All(s => !s.Sequence.Contains("*")));
    }

    #endregion

    #region Record Types Tests

    [Test]
    public void SequenceRecord_Properties_Work()
    {
        var record = new SequenceIO.SequenceRecord(
            Id: "test",
            Accession: "ACC001",
            Description: "desc",
            Sequence: "ATGC",
            Organism: "org",
            Taxonomy: "tax",
            Date: new DateTime(2024, 1, 1),
            Features: new List<SequenceIO.SequenceFeature>(),
            References: new List<SequenceIO.Reference>(),
            Metadata: new Dictionary<string, string>());

        Assert.That(record.Id, Is.EqualTo("test"));
        Assert.That(record.Sequence, Is.EqualTo("ATGC"));
    }

    [Test]
    public void BedEntry_Properties_Work()
    {
        var entry = new SequenceIO.BedEntry(
            Chromosome: "chr1",
            Start: 100,
            End: 200,
            Name: "gene1",
            Score: 500,
            Strand: '+');

        Assert.That(entry.Chromosome, Is.EqualTo("chr1"));
        Assert.That(entry.Score, Is.EqualTo(500));
    }

    [Test]
    public void GffEntry_Properties_Work()
    {
        var entry = new SequenceIO.GffEntry(
            SeqId: "chr1",
            Source: "src",
            Type: "gene",
            Start: 100,
            End: 500,
            Score: 50.5,
            Strand: '+',
            Phase: 0,
            Attributes: new Dictionary<string, string> { { "ID", "g1" } });

        Assert.That(entry.SeqId, Is.EqualTo("chr1"));
        Assert.That(entry.Score, Is.EqualTo(50.5));
    }

    [Test]
    public void SamRecord_Properties_Work()
    {
        var record = new SequenceIO.SamRecord(
            ReadName: "read1",
            Flag: 99,
            ReferenceName: "chr1",
            Position: 100,
            MappingQuality: 60,
            Cigar: "50M",
            Sequence: "ATGC",
            Quality: "IIII");

        Assert.That(record.ReadName, Is.EqualTo("read1"));
        Assert.That(record.Cigar, Is.EqualTo("50M"));
    }

    [Test]
    public void VcfRecord_Properties_Work()
    {
        var record = new SequenceIO.VcfRecord(
            Chromosome: "chr1",
            Position: 100,
            Id: "rs123",
            Reference: "A",
            Alternatives: new List<string> { "G" },
            Quality: 30.0,
            Filter: "PASS",
            Info: new Dictionary<string, string> { { "DP", "50" } },
            SampleData: new List<string> { "0/1" });

        Assert.That(record.Chromosome, Is.EqualTo("chr1"));
        Assert.That(record.Reference, Is.EqualTo("A"));
    }

    [Test]
    public void Reference_Properties_Work()
    {
        var reference = new SequenceIO.Reference(
            Number: 1,
            Authors: "Author A",
            Title: "Test Title",
            Journal: "Test Journal",
            PubMed: "12345");

        Assert.That(reference.Number, Is.EqualTo(1));
        Assert.That(reference.Authors, Is.EqualTo("Author A"));
    }

    [Test]
    public void SequenceFeature_Properties_Work()
    {
        var feature = new SequenceIO.SequenceFeature(
            Type: "gene",
            Start: 1,
            End: 100,
            Strand: '+',
            Qualifiers: new Dictionary<string, string> { { "gene", "test" } });

        Assert.That(feature.Type, Is.EqualTo("gene"));
        Assert.That(feature.Qualifiers["gene"], Is.EqualTo("test"));
    }

    #endregion
}
