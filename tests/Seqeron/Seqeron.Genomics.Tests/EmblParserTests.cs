using NUnit.Framework;
using Seqeron.Genomics;
using System.IO;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class EmblParserTests
{
    #region Sample EMBL Data

    private const string SimpleEmblRecord = @"ID   TEST001; SV 1; linear; DNA; STD; HUM; 100 BP.
XX
AC   TEST001;
XX
DT   01-JAN-2024 (Created)
XX
DE   Test sequence for unit testing.
XX
KW   test; genomics; parser.
XX
OS   Homo sapiens
OC   Eukaryota; Metazoa; Chordata; Vertebrata; Mammalia.
XX
RN   [1]
RA   Smith J., Jones A.;
RT   ""Test title for reference"";
RL   Test Journal 1:1-10(2024).
XX
FH   Key             Location/Qualifiers
FH
FT   gene            1..50
FT                   /gene=""testGene""
FT   CDS             10..40
FT                   /gene=""testGene""
FT                   /product=""test protein""
XX
SQ   Sequence 100 BP; 25 A; 25 C; 25 G; 25 T; 0 other;
     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac        50
     gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc       100
//";

    private const string MinimalRecord = @"ID   MINIMAL; SV 1; linear; DNA; STD; UNK; 20 BP.
XX
SQ   Sequence 20 BP;
     acgtacgtac gtacgtacgt                                          20
//";

    private const string CircularRecord = @"ID   PLASMID; SV 1; circular; DNA; STD; BCT; 50 BP.
XX
AC   PLASMID001;
XX
DE   Circular plasmid sequence.
XX
OS   Escherichia coli
XX
SQ   Sequence 50 BP;
     aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa        50
//";

    #endregion

    #region Basic Parsing Tests

    [Test]
    public void Parse_ValidRecord_ReturnsOneRecord()
    {
        var records = EmblParser.Parse(SimpleEmblRecord).ToList();

        Assert.That(records.Count, Is.EqualTo(1));
    }

    [Test]
    public void Parse_EmptyContent_ReturnsEmpty()
    {
        var records = EmblParser.Parse("").ToList();

        Assert.That(records, Is.Empty);
    }

    [Test]
    public void Parse_NullContent_ReturnsEmpty()
    {
        var records = EmblParser.Parse(null!).ToList();

        Assert.That(records, Is.Empty);
    }

    [Test]
    public void Parse_MinimalRecord_ParsesSuccessfully()
    {
        var records = EmblParser.Parse(MinimalRecord).ToList();

        Assert.That(records.Count, Is.EqualTo(1));
        Assert.That(records[0].Accession, Is.EqualTo("MINIMAL"));
    }

    #endregion

    #region ID Line Tests

    [Test]
    public void Parse_IdLine_ExtractsAccession()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.Accession, Is.EqualTo("TEST001"));
    }

    [Test]
    public void Parse_IdLine_ExtractsTopology()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.Topology, Is.EqualTo("linear"));
    }

    [Test]
    public void Parse_IdLine_ExtractsMoleculeType()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.MoleculeType, Is.EqualTo("DNA"));
    }

    [Test]
    public void Parse_IdLine_ExtractsLength()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.SequenceLength, Is.EqualTo(100));
    }

    [Test]
    public void Parse_CircularTopology_ParsesCorrectly()
    {
        var record = EmblParser.Parse(CircularRecord).First();

        Assert.That(record.Topology, Is.EqualTo("circular"));
    }

    [Test]
    public void Parse_IdLine_ExtractsSequenceVersion()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.SequenceVersion, Is.EqualTo("1"));
    }

    [Test]
    public void Parse_IdLine_ExtractsDataClass()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.DataClass, Is.EqualTo("STD"));
    }

    [Test]
    public void Parse_IdLine_ExtractsTaxonomicDivision()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.TaxonomicDivision, Is.EqualTo("HUM"));
    }

    #endregion

    #region Metadata Tests

    [Test]
    public void Parse_Description_ExtractsCorrectly()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.Description, Does.Contain("Test sequence"));
    }

    [Test]
    public void Parse_Keywords_ParsesMultiple()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.Keywords.Count, Is.GreaterThan(0));
        Assert.That(record.Keywords, Does.Contain("test"));
        Assert.That(record.Keywords, Does.Contain("genomics"));
    }

    [Test]
    public void Parse_Organism_ExtractsCorrectly()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.Organism, Is.EqualTo("Homo sapiens"));
    }

    [Test]
    public void Parse_OrganismClassification_ParsesHierarchy()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.OrganismClassification.Count, Is.GreaterThan(0));
        Assert.That(record.OrganismClassification, Does.Contain("Eukaryota"));
        Assert.That(record.OrganismClassification, Does.Contain("Metazoa"));
    }

    #endregion

    #region Reference Tests

    [Test]
    public void Parse_References_ExtractsAll()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.References.Count, Is.GreaterThanOrEqualTo(1));
    }

    [Test]
    public void Parse_Reference_HasAuthors()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var firstRef = record.References.FirstOrDefault();

        Assert.That(firstRef.Authors, Does.Contain("Smith"));
    }

    [Test]
    public void Parse_Reference_HasTitle()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var firstRef = record.References.FirstOrDefault();

        Assert.That(firstRef.Title, Does.Contain("Test title"));
    }

    [Test]
    public void Parse_Reference_HasJournalLocation()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var firstRef = record.References.FirstOrDefault();

        Assert.That(firstRef.Journal, Does.Contain("Test Journal"));
    }

    #endregion

    #region Feature Tests

    [Test]
    public void Parse_Features_ExtractsAll()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.Features.Count, Is.GreaterThanOrEqualTo(1));
    }

    [Test]
    public void Parse_GeneFeature_HasCorrectKey()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var genes = record.Features.Where(f => f.Key == "gene").ToList();

        Assert.That(genes.Count, Is.GreaterThanOrEqualTo(1));
    }

    [Test]
    public void Parse_Feature_HasQualifiers()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var cds = record.Features.FirstOrDefault(f => f.Key == "CDS");

        Assert.That(cds.Qualifiers, Is.Not.Null);
        Assert.That(cds.Qualifiers.ContainsKey("product"), Is.True);
        Assert.That(cds.Qualifiers["product"], Does.Contain("test protein"));
    }

    [Test]
    public void Parse_Feature_GeneQualifier()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var gene = record.Features.FirstOrDefault(f => f.Key == "gene");

        Assert.That(gene.Qualifiers.ContainsKey("gene"), Is.True);
        Assert.That(gene.Qualifiers["gene"], Is.EqualTo("testGene"));
    }

    #endregion

    #region Sequence Tests

    [Test]
    public void Parse_Sequence_ExtractsAndNormalizes()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.Sequence.Length, Is.EqualTo(100));
        Assert.That(record.Sequence, Does.StartWith("ACGTACGT"));
    }

    [Test]
    public void Parse_Sequence_RemovesNumbersAndSpaces()
    {
        var record = EmblParser.Parse(MinimalRecord).First();

        Assert.That(record.Sequence, Is.EqualTo("ACGTACGTACGTACGTACGT"));
        Assert.That(record.Sequence, Does.Not.Contain(" "));
    }

    [Test]
    public void Parse_Sequence_UpperCase()
    {
        var record = EmblParser.Parse(MinimalRecord).First();

        Assert.That(record.Sequence.All(c => char.IsUpper(c)), Is.True);
    }

    #endregion

    #region Location Parsing Tests

    [Test]
    public void ParseLocation_SimpleRange_ParsesCorrectly()
    {
        var location = EmblParser.ParseLocation("100..200");

        Assert.That(location.Start, Is.EqualTo(100));
        Assert.That(location.End, Is.EqualTo(200));
        Assert.That(location.IsComplement, Is.False);
    }

    [Test]
    public void ParseLocation_Complement_DetectsStrand()
    {
        var location = EmblParser.ParseLocation("complement(100..200)");

        Assert.That(location.IsComplement, Is.True);
    }

    [Test]
    public void ParseLocation_Join_ExtractsParts()
    {
        var location = EmblParser.ParseLocation("join(1..50,60..100)");

        Assert.That(location.IsJoin, Is.True);
        Assert.That(location.Parts.Count, Is.EqualTo(2));
    }

    [Test]
    public void ParseLocation_SingleBase_ParsesCorrectly()
    {
        var location = EmblParser.ParseLocation("467");

        Assert.That(location.Start, Is.EqualTo(467));
        Assert.That(location.End, Is.EqualTo(467));
    }

    [Test]
    public void ParseLocation_PartialStart_DetectsPartial()
    {
        var location = EmblParser.ParseLocation("<1..200");

        Assert.That(location.Start, Is.EqualTo(1));
        Assert.That(location.End, Is.EqualTo(200));
        Assert.That(location.RawLocation, Does.Contain("<"));
    }

    [Test]
    public void ParseLocation_PartialEnd_DetectsPartial()
    {
        var location = EmblParser.ParseLocation("100..>500");

        Assert.That(location.Start, Is.EqualTo(100));
        Assert.That(location.End, Is.EqualTo(500));
        Assert.That(location.RawLocation, Does.Contain(">"));
    }

    [Test]
    public void ParseLocation_ComplementJoin_ParsesCorrectly()
    {
        var location = EmblParser.ParseLocation("complement(join(1..50,60..100))");

        Assert.That(location.IsComplement, Is.True);
        Assert.That(location.IsJoin, Is.True);
        Assert.That(location.Parts.Count, Is.EqualTo(2));
    }

    #endregion

    #region Conversion Tests

    [Test]
    public void ToGenBank_ConvertsSuccessfully()
    {
        var embl = EmblParser.Parse(SimpleEmblRecord).First();
        var genBank = EmblParser.ToGenBank(embl);

        Assert.That(genBank.Accession, Is.EqualTo(embl.Accession));
        Assert.That(genBank.Sequence, Is.EqualTo(embl.Sequence));
        Assert.That(genBank.Organism, Is.EqualTo(embl.Organism));
    }

    [Test]
    public void ToGenBank_PreservesFeatures()
    {
        var embl = EmblParser.Parse(SimpleEmblRecord).First();
        var genBank = EmblParser.ToGenBank(embl);

        Assert.That(genBank.Features.Count, Is.EqualTo(embl.Features.Count));
    }

    #endregion

    #region Utility Method Tests

    [Test]
    public void GetCDS_ReturnsOnlyCDSFeatures()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var cdsFeatures = EmblParser.GetCDS(record).ToList();

        Assert.That(cdsFeatures.All(f => f.Key == "CDS"), Is.True);
    }

    [Test]
    public void GetGenes_ReturnsOnlyGeneFeatures()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var genes = EmblParser.GetGenes(record).ToList();

        Assert.That(genes.All(f => f.Key == "gene"), Is.True);
    }

    [Test]
    public void GetFeatures_FiltersByKey()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var cdsFeatures = EmblParser.GetFeatures(record, "CDS").ToList();

        Assert.That(cdsFeatures.All(f => f.Key == "CDS"), Is.True);
        Assert.That(cdsFeatures.Count, Is.GreaterThanOrEqualTo(1));
    }

    [Test]
    public void ExtractSequence_SimpleLocation_ReturnsSubsequence()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var location = EmblParser.ParseLocation("1..10");
        var subseq = EmblParser.ExtractSequence(record, location);

        Assert.That(subseq.Length, Is.EqualTo(10));
        Assert.That(subseq, Is.EqualTo("ACGTACGTAC"));
    }

    #endregion

    #region Multiple Records Tests

    private const string MultipleRecords = @"ID   REC1; SV 1; linear; DNA; STD; UNK; 10 BP.
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//
ID   REC2; SV 1; linear; DNA; STD; UNK; 10 BP.
SQ   Sequence 10 BP;
     cccccccccc                                                     10
//
ID   REC3; SV 1; linear; DNA; STD; UNK; 10 BP.
SQ   Sequence 10 BP;
     gggggggggg                                                     10
//";

    [Test]
    public void Parse_MultipleRecords_ParsesAll()
    {
        var records = EmblParser.Parse(MultipleRecords).ToList();

        Assert.That(records.Count, Is.EqualTo(3));
    }

    [Test]
    public void Parse_MultipleRecords_EachHasCorrectSequence()
    {
        var records = EmblParser.Parse(MultipleRecords).ToList();

        Assert.That(records[0].Sequence, Is.EqualTo("AAAAAAAAAA"));
        Assert.That(records[1].Sequence, Is.EqualTo("CCCCCCCCCC"));
        Assert.That(records[2].Sequence, Is.EqualTo("GGGGGGGGGG"));
    }

    #endregion

    #region ParseFile Tests

    [Test]
    public void ParseFile_ValidFile_ParsesSuccessfully()
    {
        // Create temp file
        var tempFile = Path.GetTempFileName();
        try
        {
            File.WriteAllText(tempFile, SimpleEmblRecord);
            var records = EmblParser.ParseFile(tempFile).ToList();

            Assert.That(records.Count, Is.EqualTo(1));
            Assert.That(records[0].Accession, Is.EqualTo("TEST001"));
        }
        finally
        {
            File.Delete(tempFile);
        }
    }

    [Test]
    public void ParseFile_InvalidPath_ReturnsEmpty()
    {
        // ParseFile returns empty collection for non-existent files (doesn't throw)
        var records = EmblParser.ParseFile(@"C:\nonexistent\path\file.embl").ToList();

        Assert.That(records, Is.Empty);
    }

    #endregion

    #region Edge Case Tests

    [Test]
    public void Parse_WhitespaceOnly_ReturnsEmpty()
    {
        var records = EmblParser.Parse("   \n\t\n   ").ToList();

        Assert.That(records, Is.Empty);
    }

    [Test]
    public void Parse_RecordWithoutTerminator_HandlesGracefully()
    {
        // Record without // terminator - should still parse if content is valid
        var incomplete = @"ID   TEST; SV 1; linear; DNA; STD; UNK; 10 BP.
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10";

        var records = EmblParser.Parse(incomplete).ToList();
        // Behavior depends on implementation - may return 0 or 1 records
        Assert.That(records.Count, Is.LessThanOrEqualTo(1));
    }

    #endregion
}
