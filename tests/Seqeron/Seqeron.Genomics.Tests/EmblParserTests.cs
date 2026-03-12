using NUnit.Framework;
using Seqeron.Genomics;
using System.IO;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class EmblParserTests
{
    #region Sample EMBL Data

    private const string SimpleEmblRecord = @"ID   TEST001; SV 1; linear; genomic DNA; STD; HUM; 100 BP.
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

    private const string MinimalRecord = @"ID   MINIMAL; SV 1; linear; genomic DNA; STD; UNC; 20 BP.
XX
SQ   Sequence 20 BP;
     acgtacgtac gtacgtacgt                                          20
//";

    private const string CircularRecord = @"ID   PLASMID; SV 1; circular; genomic DNA; STD; PRO; 50 BP.
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

        Assert.That(record.MoleculeType, Is.EqualTo("genomic DNA"));
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

        Assert.That(record.Description, Is.EqualTo("Test sequence for unit testing."));
    }

    [Test]
    public void Parse_Keywords_ParsesMultiple()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.Keywords.Count, Is.EqualTo(3));
        Assert.That(record.Keywords, Is.EqualTo(new[] { "test", "genomics", "parser" }));
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

        Assert.That(record.OrganismClassification.Count, Is.EqualTo(5));
        Assert.That(record.OrganismClassification, Is.EqualTo(new[] { "Eukaryota", "Metazoa", "Chordata", "Vertebrata", "Mammalia" }));
    }

    #endregion

    #region Reference Tests

    [Test]
    public void Parse_References_ExtractsAll()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.References.Count, Is.EqualTo(1));
    }

    [Test]
    public void Parse_Reference_HasAuthors()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var firstRef = record.References[0];

        Assert.That(firstRef.Authors, Is.EqualTo("Smith J., Jones A."));
    }

    [Test]
    public void Parse_Reference_HasTitle()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var firstRef = record.References[0];

        Assert.That(firstRef.Title, Is.EqualTo("Test title for reference"));
    }

    [Test]
    public void Parse_Reference_HasJournalLocation()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var firstRef = record.References[0];

        Assert.That(firstRef.Journal, Is.EqualTo("Test Journal 1:1-10(2024)."));
    }

    #endregion

    #region Feature Tests

    [Test]
    public void Parse_Features_ExtractsAll()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();

        Assert.That(record.Features.Count, Is.EqualTo(2));
        Assert.That(record.Features[0].Key, Is.EqualTo("gene"));
        Assert.That(record.Features[1].Key, Is.EqualTo("CDS"));
    }

    [Test]
    public void Parse_GeneFeature_HasCorrectKey()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var genes = record.Features.Where(f => f.Key == "gene").ToList();

        Assert.That(genes.Count, Is.EqualTo(1));
        Assert.That(genes[0].Location.Start, Is.EqualTo(1));
        Assert.That(genes[0].Location.End, Is.EqualTo(50));
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
        Assert.That(location.Start, Is.EqualTo(100));
        Assert.That(location.End, Is.EqualTo(200));
        Assert.That(location.IsJoin, Is.False);
    }

    [Test]
    public void ParseLocation_Join_ExtractsParts()
    {
        var location = EmblParser.ParseLocation("join(1..50,60..100)");

        Assert.That(location.IsJoin, Is.True);
        Assert.That(location.IsComplement, Is.False);
        Assert.That(location.Parts.Count, Is.EqualTo(2));
        Assert.That(location.Parts[0], Is.EqualTo((1, 50)));
        Assert.That(location.Parts[1], Is.EqualTo((60, 100)));
        Assert.That(location.Start, Is.EqualTo(1));
        Assert.That(location.End, Is.EqualTo(100));
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
        Assert.That(location.Is5PrimePartial, Is.True);
        Assert.That(location.Is3PrimePartial, Is.False);
        Assert.That(location.RawLocation, Is.EqualTo("<1..200"));
    }

    [Test]
    public void ParseLocation_PartialEnd_DetectsPartial()
    {
        var location = EmblParser.ParseLocation("100..>500");

        Assert.That(location.Start, Is.EqualTo(100));
        Assert.That(location.End, Is.EqualTo(500));
        Assert.That(location.Is3PrimePartial, Is.True);
        Assert.That(location.Is5PrimePartial, Is.False);
        Assert.That(location.RawLocation, Is.EqualTo("100..>500"));
    }

    [Test]
    public void ParseLocation_ComplementJoin_ParsesCorrectly()
    {
        var location = EmblParser.ParseLocation("complement(join(1..50,60..100))");

        Assert.That(location.IsComplement, Is.True);
        Assert.That(location.IsJoin, Is.True);
        Assert.That(location.Parts.Count, Is.EqualTo(2));
        Assert.That(location.Parts[0], Is.EqualTo((1, 50)));
        Assert.That(location.Parts[1], Is.EqualTo((60, 100)));
        Assert.That(location.Start, Is.EqualTo(1));
        Assert.That(location.End, Is.EqualTo(100));
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

        Assert.That(cdsFeatures.Count, Is.EqualTo(1));
        Assert.That(cdsFeatures[0].Key, Is.EqualTo("CDS"));
        Assert.That(cdsFeatures[0].Location.Start, Is.EqualTo(10));
        Assert.That(cdsFeatures[0].Location.End, Is.EqualTo(40));
    }

    [Test]
    public void GetGenes_ReturnsOnlyGeneFeatures()
    {
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var genes = EmblParser.GetGenes(record).ToList();

        Assert.That(genes.Count, Is.EqualTo(1));
        Assert.That(genes[0].Key, Is.EqualTo("gene"));
        Assert.That(genes[0].Qualifiers["gene"], Is.EqualTo("testGene"));
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

    private const string MultipleRecords = @"ID   REC1; SV 1; linear; genomic DNA; STD; UNC; 10 BP.
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//
ID   REC2; SV 1; linear; genomic DNA; STD; UNC; 10 BP.
SQ   Sequence 10 BP;
     cccccccccc                                                     10
//
ID   REC3; SV 1; linear; genomic DNA; STD; UNC; 10 BP.
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
        // Record without // terminator - parser splits by \n// so unterminated content
        // remains as a single chunk and is parsed if it starts with ID.
        var incomplete = @"ID   TEST; SV 1; linear; genomic DNA; STD; UNC; 10 BP.
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10";

        var records = EmblParser.Parse(incomplete).ToList();

        Assert.That(records.Count, Is.EqualTo(1));
        Assert.That(records[0].Accession, Is.EqualTo("TEST"));
        Assert.That(records[0].Sequence, Is.EqualTo("AAAAAAAAAA"));
    }

    #endregion

    #region Location Parsing — Order

    [Test]
    public void ParseLocation_Order_ParsesCorrectly()
    {
        var location = EmblParser.ParseLocation("order(100..200,300..400)");

        Assert.That(location.IsOrder, Is.True);
        Assert.That(location.IsJoin, Is.False);
        Assert.That(location.IsComplement, Is.False);
        Assert.That(location.Parts.Count, Is.EqualTo(2));
        Assert.That(location.Parts[0], Is.EqualTo((100, 200)));
        Assert.That(location.Parts[1], Is.EqualTo((300, 400)));
        Assert.That(location.Start, Is.EqualTo(100));
        Assert.That(location.End, Is.EqualTo(400));
    }

    #endregion

    #region Reference — DOI and PubMed

    [Test]
    public void Parse_Reference_HasDOI()
    {
        var embl = @"ID   DOI001; SV 1; linear; mRNA; STD; PLN; 10 BP.
XX
RN   [1]
RX   DOI; 10.1007/BF00039495.
RA   Smith J.;
RT   ""Title"";
RL   Journal 1:1-10(2024).
XX
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.References[0].CrossReference, Does.Contain("10.1007/BF00039495"));
    }

    [Test]
    public void Parse_Reference_HasPubMed()
    {
        var embl = @"ID   PM001; SV 1; linear; mRNA; STD; PLN; 10 BP.
XX
RN   [1]
RX   PUBMED; 1907511.
RA   Smith J.;
RT   ""Title"";
RL   Journal 1:1-10(2024).
XX
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.References[0].CrossReference, Does.Contain("1907511"));
    }

    #endregion

    #region AdditionalFields — OG, DR, CC, DT

    [Test]
    public void Parse_Organelle_ExtractsCorrectly()
    {
        // OG line: organelle information — stored in AdditionalFields["OG"].
        var embl = @"ID   ORG001; SV 1; linear; genomic DNA; STD; PLN; 10 BP.
XX
OG   Mitochondrion
XX
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.AdditionalFields.ContainsKey("OG"), Is.True);
        Assert.That(record.AdditionalFields["OG"], Is.EqualTo("Mitochondrion"));
    }

    [Test]
    public void Parse_DatabaseCrossReference_DR()
    {
        // DR line: database cross-reference — stored in AdditionalFields["DR"].
        var embl = @"ID   DR001; SV 1; linear; genomic DNA; STD; HUM; 10 BP.
XX
DR   UniProtKB/Swiss-Prot; P26204; AMYG_TRIRP.
XX
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.AdditionalFields.ContainsKey("DR"), Is.True);
        Assert.That(record.AdditionalFields["DR"], Does.Contain("UniProtKB/Swiss-Prot"));
        Assert.That(record.AdditionalFields["DR"], Does.Contain("P26204"));
    }

    [Test]
    public void Parse_Comments_CC()
    {
        // CC line: free-text comments — stored in AdditionalFields["CC"].
        var embl = @"ID   CC001; SV 1; linear; genomic DNA; STD; HUM; 10 BP.
XX
CC   This is a comment about the sequence.
XX
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.AdditionalFields.ContainsKey("CC"), Is.True);
        Assert.That(record.AdditionalFields["CC"], Is.EqualTo("This is a comment about the sequence."));
    }

    [Test]
    public void Parse_DateLines_DT()
    {
        // DT line: date information — stored in AdditionalFields["DT"].
        var embl = @"ID   DT001; SV 1; linear; genomic DNA; STD; HUM; 10 BP.
XX
DT   01-JAN-2024 (Created)
XX
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.AdditionalFields.ContainsKey("DT"), Is.True);
        Assert.That(record.AdditionalFields["DT"], Does.Contain("01-JAN-2024"));
    }

    #endregion

    #region Multi-line Continuation and Edge Cases

    [Test]
    public void Parse_MultiLineContinuation_DE()
    {
        // DE lines can span multiple lines; content should be joined.
        var embl = @"ID   ML001; SV 1; linear; genomic DNA; STD; HUM; 10 BP.
XX
DE   This is a long description that spans
DE   multiple lines in the EMBL file.
XX
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.Description, Is.EqualTo("This is a long description that spans multiple lines in the EMBL file."));
    }

    [Test]
    public void Parse_EmptyKeywords_ReturnsEmpty()
    {
        // "KW   ." indicates no keywords.
        var embl = @"ID   EK001; SV 1; linear; genomic DNA; STD; HUM; 10 BP.
XX
KW   .
XX
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.Keywords, Is.Empty);
    }

    [Test]
    public void Parse_SecondaryAccessions_PrimaryExtracted()
    {
        // AC line may contain primary and secondary accessions separated by semicolons.
        // Parser extracts primary accession; secondaries are not stored individually.
        var embl = @"ID   PRIM01; SV 1; linear; genomic DNA; STD; HUM; 10 BP.
XX
AC   PRIM01; SEC001; SEC002;
XX
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.Accession, Is.EqualTo("PRIM01"));
    }

    #endregion

    #region ExtractSequence — Complement and Join

    [Test]
    public void ExtractSequence_ComplementLocation_ReturnsReverseComplement()
    {
        // SimpleEmblRecord sequence starts with "ACGTACGTAC..." (100 chars).
        // complement(1..10) → extract "ACGTACGTAC", then reverse complement → "GTACGTACGT".
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var location = EmblParser.ParseLocation("complement(1..10)");
        var subseq = EmblParser.ExtractSequence(record, location);

        Assert.That(subseq, Is.EqualTo("GTACGTACGT"));
    }

    [Test]
    public void ExtractSequence_JoinLocation_ReturnsConcatenatedRegions()
    {
        // SimpleEmblRecord: positions 1..5 = "ACGTA", positions 51..55 = "GCGCG".
        // join(1..5,51..55) → "ACGTAGCGCG".
        var record = EmblParser.Parse(SimpleEmblRecord).First();
        var location = EmblParser.ParseLocation("join(1..5,51..55)");
        var subseq = EmblParser.ExtractSequence(record, location);

        Assert.That(subseq, Is.EqualTo("ACGTAGCGCG"));
    }

    #endregion

    #region Spec Compliance Tests — INSDC Feature Table v11.3 & EBI EMBL User Manual Release 143

    [TestCase("genomic DNA")]
    [TestCase("genomic RNA")]
    [TestCase("mRNA")]
    [TestCase("tRNA")]
    [TestCase("rRNA")]
    [TestCase("other RNA")]
    [TestCase("other DNA")]
    [TestCase("transcribed RNA")]
    [TestCase("viral cRNA")]
    [TestCase("unassigned DNA")]
    [TestCase("unassigned RNA")]
    public void Parse_IdLine_AllInsdcMolTypes(string molType)
    {
        var embl = $@"ID   MOL001; SV 1; linear; {molType}; STD; HUM; 10 BP.
XX
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.MoleculeType, Is.EqualTo(molType));
    }

    [Test]
    public void Parse_IdLine_BareDnaNotRecognisedAsMolType()
    {
        // Bare "DNA" is not in the INSDC mol_type vocabulary; field should be empty.
        var embl = @"ID   BARE; SV 1; linear; DNA; STD; HUM; 10 BP.
XX
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.MoleculeType, Is.Empty);
    }

    [TestCase("CON")]
    [TestCase("PAT")]
    [TestCase("EST")]
    [TestCase("GSS")]
    [TestCase("HTC")]
    [TestCase("HTG")]
    [TestCase("WGS")]
    [TestCase("TSA")]
    [TestCase("STS")]
    [TestCase("STD")]
    public void Parse_IdLine_AllDataClasses(string dataClass)
    {
        var embl = $@"ID   DC001; SV 1; linear; genomic DNA; {dataClass}; HUM; 10 BP.
XX
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.DataClass, Is.EqualTo(dataClass));
    }

    [TestCase("PHG")]
    [TestCase("ENV")]
    [TestCase("FUN")]
    [TestCase("HUM")]
    [TestCase("INV")]
    [TestCase("MAM")]
    [TestCase("VRT")]
    [TestCase("MUS")]
    [TestCase("PLN")]
    [TestCase("PRO")]
    [TestCase("ROD")]
    [TestCase("SYN")]
    [TestCase("TGN")]
    [TestCase("UNC")]
    [TestCase("VRL")]
    public void Parse_IdLine_AllTaxonomicDivisions(string division)
    {
        var embl = $@"ID   DIV001; SV 1; linear; genomic DNA; STD; {division}; 10 BP.
XX
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.TaxonomicDivision, Is.EqualTo(division));
    }

    [Test]
    public void Parse_IdLine_InvalidDivisionNotRecognised()
    {
        // "UNK" is not a valid EMBL division code.
        var embl = @"ID   INV001; SV 1; linear; genomic DNA; STD; UNK; 10 BP.
XX
SQ   Sequence 10 BP;
     aaaaaaaaaa                                                     10
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.TaxonomicDivision, Is.Empty);
    }

    [Test]
    public void Parse_Qualifier_SlashInValue_NotTruncated()
    {
        // Regression test: /db_xref values with "/" must not be truncated.
        // INSDC Feature Table v11.3: qualifier values are delimited by quotes, not by "/".
        var embl = @"ID   QUAL001; SV 1; linear; mRNA; STD; PLN; 20 BP.
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..20
FT                   /organism=""Trifolium repens""
FT                   /db_xref=""UniProtKB/Swiss-Prot:P26204""
XX
SQ   Sequence 20 BP;
     acgtacgtac gtacgtacgt                                          20
//";
        var record = EmblParser.Parse(embl).First();
        var source = record.Features.First(f => f.Key == "source");

        Assert.That(source.Qualifiers["db_xref"], Is.EqualTo("UniProtKB/Swiss-Prot:P26204"));
    }

    [Test]
    public void Parse_Qualifier_MultipleSlashesInValue()
    {
        // Ensure multiple slashes in a quoted value are preserved.
        var embl = @"ID   QUAL002; SV 1; linear; mRNA; STD; PLN; 20 BP.
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..20
FT                   /note=""path/to/some/resource""
XX
SQ   Sequence 20 BP;
     acgtacgtac gtacgtacgt                                          20
//";
        var record = EmblParser.Parse(embl).First();
        var source = record.Features.First(f => f.Key == "source");

        Assert.That(source.Qualifiers["note"], Is.EqualTo("path/to/some/resource"));
    }

    [Test]
    public void Parse_Reference_CapturesPositions()
    {
        // RP line: "Reference Positions" — per EBI User Manual Section 3.4.10.
        var embl = @"ID   REF001; SV 1; linear; mRNA; STD; PLN; 20 BP.
XX
RN   [1]
RP   1-1859
RA   Oxtoby E., Dunn M.A.;
RT   ""Nucleotide sequence"";
RL   Plant Mol. Biol. 17(2):209-219(1991).
XX
SQ   Sequence 20 BP;
     acgtacgtac gtacgtacgt                                          20
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.References[0].Positions, Is.EqualTo("1-1859"));
    }

    [Test]
    public void Parse_Reference_CapturesGroup()
    {
        // RG line: "Reference Group" — per EBI User Manual Section 3.4.10.
        var embl = @"ID   REF002; SV 1; linear; mRNA; STD; PLN; 20 BP.
XX
RN   [1]
RP   1-20
RG   The Genome Consortium
RA   Smith J.;
RT   ""Title"";
RL   Journal 1:1-10(2024).
XX
SQ   Sequence 20 BP;
     acgtacgtac gtacgtacgt                                          20
//";
        var record = EmblParser.Parse(embl).First();

        Assert.That(record.References[0].Group, Is.EqualTo("The Genome Consortium"));
    }

    [Test]
    public void Parse_EbiReferenceRecord_FullBlock()
    {
        // Complete reference block from EBI User Manual example (X56734).
        var embl = @"ID   X56734; SV 1; linear; mRNA; STD; PLN; 30 BP.
XX
AC   X56734; S46826;
XX
RN   [5]
RP   1-30
RX   DOI; 10.1007/BF00039495.
RX   PUBMED; 1907511.
RA   Oxtoby E., Dunn M.A., Pancoro A., Hughes M.A.;
RT   ""Nucleotide and derived amino acid sequence"";
RL   Plant Mol. Biol. 17(2):209-219(1991).
XX
SQ   Sequence 30 BP;
     acgtacgtac gtacgtacgt acgtacgtac                               30
//";
        var record = EmblParser.Parse(embl).First();
        var r = record.References[0];

        Assert.That(r.Number, Is.EqualTo(5));
        Assert.That(r.Positions, Is.EqualTo("1-30"));
        Assert.That(r.Authors, Does.Contain("Oxtoby"));
        Assert.That(r.Title, Does.Contain("Nucleotide"));
        Assert.That(r.Journal, Does.Contain("Plant Mol. Biol."));
        Assert.That(r.CrossReference, Does.Contain("1907511"));
    }

    [Test]
    public void ParseLocation_SiteBetween_ParsesRange()
    {
        // INSDC Feature Table: 123^124 means a site between bases 123 and 124.
        // The parser captures Start=123, End=124 (the flanking positions).
        var location = EmblParser.ParseLocation("123^124");

        Assert.That(location.Start, Is.EqualTo(123));
        Assert.That(location.End, Is.EqualTo(124));
        Assert.That(location.RawLocation, Is.EqualTo("123^124"));
    }

    #endregion
}
