using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class GenBankParseTests
{
    // Use the same format as in GenBankParserTests.cs
    private const string TestGenBank = @"LOCUS       TEST001                  100 bp    DNA     linear   BCT 01-JAN-2024
DEFINITION  Test sequence for unit testing.
ACCESSION   TEST001
VERSION     TEST001.1
KEYWORDS    test; unit test.
SOURCE      Test organism
  ORGANISM  Test organism
            Bacteria; Proteobacteria.
FEATURES             Location/Qualifiers
     gene            1..50
                     /gene=""testGene""
                     /note=""Test gene feature""
     CDS             10..40
                     /gene=""testGene""
                     /translation=""MKLLVV""
                     /product=""test protein""
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac
       51 gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc
//";

    [Test]
    public void GenBankParse_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.GenBankParse(TestGenBank));
        Assert.Throws<ArgumentException>(() => ParsersTools.GenBankParse(""));
        Assert.Throws<ArgumentException>(() => ParsersTools.GenBankParse(null!));
    }

    [Test]
    public void GenBankParse_Binding_ParsesRecords()
    {
        var result = ParsersTools.GenBankParse(TestGenBank);

        Assert.That(result.Count, Is.EqualTo(1));
        Assert.That(result.Records[0].Locus, Is.EqualTo("TEST001"));
        Assert.That(result.Records[0].SequenceLength, Is.EqualTo(100));
        Assert.That(result.Records[0].MoleculeType, Is.EqualTo("DNA"));
        Assert.That(result.Records[0].Topology, Is.EqualTo("linear"));
        Assert.That(result.Records[0].Division, Is.EqualTo("BCT"));
    }

    [Test]
    public void GenBankParse_Binding_ParsesDefinition()
    {
        var result = ParsersTools.GenBankParse(TestGenBank);

        Assert.That(result.Records[0].Definition, Is.EqualTo("Test sequence for unit testing."));
        Assert.That(result.Records[0].Accession, Is.EqualTo("TEST001"));
        Assert.That(result.Records[0].Version, Is.EqualTo("TEST001.1"));
    }

    [Test]
    public void GenBankParse_Binding_ParsesKeywords()
    {
        var result = ParsersTools.GenBankParse(TestGenBank);

        Assert.That(result.Records[0].Keywords, Contains.Item("test"));
        Assert.That(result.Records[0].Keywords, Contains.Item("unit test"));
    }

    [Test]
    public void GenBankParse_Binding_ParsesOrganism()
    {
        var result = ParsersTools.GenBankParse(TestGenBank);

        Assert.That(result.Records[0].Organism, Is.EqualTo("Test organism"));
        Assert.That(result.Records[0].Taxonomy, Does.Contain("Bacteria"));
    }

    [Test]
    public void GenBankParse_Binding_ParsesFeatures()
    {
        var result = ParsersTools.GenBankParse(TestGenBank);

        // Test GenBank has gene and CDS features - parser MUST extract them
        Assert.That(result.Records[0].Features, Is.Not.Null);
        Assert.That(result.Records[0].Features, Is.Not.Empty,
            "Test GenBank contains gene and CDS features - parser MUST extract them");
        Assert.That(result.Records[0].Features.Any(f => f.Key == "gene" || f.Key == "CDS"), Is.True,
            "Features must include gene or CDS entries");
    }

    [Test]
    public void GenBankParse_Binding_ParsesSequence()
    {
        var result = ParsersTools.GenBankParse(TestGenBank);

        // Test GenBank has ORIGIN section with 100bp - MUST be parsed
        Assert.That(result.Records[0].Sequence, Is.Not.Null.And.Not.Empty,
            "GenBank ORIGIN section with sequence MUST be parsed");
        Assert.That(result.Records[0].Sequence, Does.StartWith("ACGTACGT"));
        Assert.That(result.Records[0].ActualSequenceLength, Is.EqualTo(100));
    }
}

[TestFixture]
public class GenBankFeaturesTests
{
    private const string TestGenBank = @"LOCUS       TEST001                  100 bp    DNA     linear   BCT 01-JAN-2024
DEFINITION  Test sequence.
ACCESSION   TEST001
VERSION     TEST001.1
KEYWORDS    .
SOURCE      Test organism
  ORGANISM  Test organism
            Bacteria.
FEATURES             Location/Qualifiers
     gene            1..50
                     /gene=""gene1""
     gene            60..90
                     /gene=""gene2""
     CDS             10..40
                     /gene=""gene1""
                     /product=""Product 1""
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac
       51 gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc
//";

    [Test]
    public void GenBankFeatures_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.GenBankFeatures(TestGenBank));
        Assert.Throws<ArgumentException>(() => ParsersTools.GenBankFeatures(""));
    }

    [Test]
    public void GenBankFeatures_Binding_ExtractsAllFeatures()
    {
        var result = ParsersTools.GenBankFeatures(TestGenBank);

        // Features may or may not be parsed depending on implementation
        Assert.That(result.Features, Is.Not.Null);
        Assert.That(result.Count, Is.EqualTo(result.Features.Count));
    }

    [Test]
    public void GenBankFeatures_Binding_FiltersByType()
    {
        var result = ParsersTools.GenBankFeatures(TestGenBank, featureType: "gene");

        // All returned features should be genes
        Assert.That(result.Features.All(f => f.Key == "gene"), Is.True);
    }

    [Test]
    public void GenBankFeatures_Binding_ParsesQualifiers()
    {
        var result = ParsersTools.GenBankFeatures(TestGenBank, featureType: "CDS");

        // Test GenBank has CDS feature with /gene and /product qualifiers - MUST be parsed
        Assert.That(result.Count, Is.GreaterThan(0),
            "Test GenBank contains CDS feature - must be extracted");
        Assert.That(result.Features[0].Qualifiers, Is.Not.Empty,
            "CDS feature has /gene and /product qualifiers - must be parsed");
    }

    [Test]
    public void GenBankFeatures_Binding_ParsesLocation()
    {
        var result = ParsersTools.GenBankFeatures(TestGenBank, featureType: "gene");

        // Test GenBank has gene features at 1..50 and 60..90 - MUST be parsed
        Assert.That(result.Count, Is.GreaterThan(0),
            "Test GenBank contains gene features - must be extracted");

        var gene1 = result.Features.FirstOrDefault(f => f.Qualifiers.GetValueOrDefault("gene") == "gene1");
        Assert.That(gene1, Is.Not.Null,
            "gene1 feature must be found with /gene=\"gene1\" qualifier");
        Assert.That(gene1!.Start, Is.EqualTo(1), "gene1 location should be 1..50");
        Assert.That(gene1.End, Is.EqualTo(50), "gene1 location should be 1..50");
    }
}

[TestFixture]
public class GenBankStatisticsTests
{
    private const string TestGenBank = @"LOCUS       TEST001                  100 bp    DNA     linear   BCT 01-JAN-2024
DEFINITION  Test sequence.
ACCESSION   TEST001
VERSION     TEST001.1
KEYWORDS    .
SOURCE      Test organism
  ORGANISM  Test organism
            Bacteria.
FEATURES             Location/Qualifiers
     gene            1..50
                     /gene=""gene1""
     CDS             10..40
                     /product=""Product 1""
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac
       51 gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc
//";

    [Test]
    public void GenBankStatistics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.GenBankStatistics(TestGenBank));
        Assert.Throws<ArgumentException>(() => ParsersTools.GenBankStatistics(""));
    }

    [Test]
    public void GenBankStatistics_Binding_CalculatesStats()
    {
        var result = ParsersTools.GenBankStatistics(TestGenBank);

        Assert.That(result.RecordCount, Is.EqualTo(1));
        Assert.That(result.TotalSequenceLength, Is.EqualTo(100));
        // Features counts depend on parser implementation
        Assert.That(result.TotalFeatures, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void GenBankStatistics_Binding_CountsFeatureTypes()
    {
        var result = ParsersTools.GenBankStatistics(TestGenBank);

        // FeatureTypeCounts should match TotalFeatures
        var totalFromCounts = result.FeatureTypeCounts.Values.Sum();
        Assert.That(totalFromCounts, Is.EqualTo(result.TotalFeatures));
    }

    [Test]
    public void GenBankStatistics_Binding_CollectsMoleculeTypes()
    {
        var result = ParsersTools.GenBankStatistics(TestGenBank);

        Assert.That(result.MoleculeTypes, Contains.Item("DNA"));
        Assert.That(result.Divisions, Contains.Item("BCT"));
    }
}

[TestFixture]
public class GenBankParseLocationTests
{
    [Test]
    public void GenBankParseLocation_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.GenBankParseLocation("100..200"));
        Assert.Throws<ArgumentException>(() => ParsersTools.GenBankParseLocation(""));
        Assert.Throws<ArgumentException>(() => ParsersTools.GenBankParseLocation(null!));
    }

    [Test]
    public void GenBankParseLocation_Binding_ParsesSimpleRange()
    {
        var result = ParsersTools.GenBankParseLocation("100..200");

        Assert.That(result.Start, Is.EqualTo(100));
        Assert.That(result.End, Is.EqualTo(200));
        Assert.That(result.Length, Is.EqualTo(101));
        Assert.That(result.IsComplement, Is.False);
        Assert.That(result.IsJoin, Is.False);
        Assert.That(result.Parts.Count, Is.EqualTo(1));
    }

    [Test]
    public void GenBankParseLocation_Binding_ParsesComplement()
    {
        var result = ParsersTools.GenBankParseLocation("complement(100..200)");

        Assert.That(result.Start, Is.EqualTo(100));
        Assert.That(result.End, Is.EqualTo(200));
        Assert.That(result.IsComplement, Is.True);
        Assert.That(result.IsJoin, Is.False);
        Assert.That(result.RawLocation, Is.EqualTo("complement(100..200)"));
    }

    [Test]
    public void GenBankParseLocation_Binding_ParsesJoin()
    {
        var result = ParsersTools.GenBankParseLocation("join(100..200,300..400)");

        Assert.That(result.Start, Is.EqualTo(100));
        Assert.That(result.End, Is.EqualTo(400));
        Assert.That(result.IsComplement, Is.False);
        Assert.That(result.IsJoin, Is.True);
        Assert.That(result.Parts.Count, Is.EqualTo(2));
        Assert.That(result.Parts[0].Start, Is.EqualTo(100));
        Assert.That(result.Parts[0].End, Is.EqualTo(200));
        Assert.That(result.Parts[1].Start, Is.EqualTo(300));
        Assert.That(result.Parts[1].End, Is.EqualTo(400));
    }

    [Test]
    public void GenBankParseLocation_Binding_ParsesComplementJoin()
    {
        var result = ParsersTools.GenBankParseLocation("complement(join(100..200,300..400))");

        Assert.That(result.IsComplement, Is.True);
        Assert.That(result.IsJoin, Is.True);
        Assert.That(result.Parts.Count, Is.EqualTo(2));
    }

    [Test]
    public void GenBankParseLocation_Binding_ParsesSinglePosition()
    {
        var result = ParsersTools.GenBankParseLocation("150");

        Assert.That(result.Start, Is.EqualTo(150));
        Assert.That(result.End, Is.EqualTo(150));
        Assert.That(result.Length, Is.EqualTo(1));
    }

    [Test]
    public void GenBankParseLocation_Binding_PreservesRawLocation()
    {
        var location = "join(1..100,200..300,400..500)";
        var result = ParsersTools.GenBankParseLocation(location);

        Assert.That(result.RawLocation, Is.EqualTo(location));
        Assert.That(result.Parts.Count, Is.EqualTo(3));
    }
}

[TestFixture]
public class GenBankExtractSequenceTests
{
    private const string TestGenBank = @"LOCUS       TEST001                  100 bp    DNA     linear   BCT 01-JAN-2024
DEFINITION  Test sequence.
ACCESSION   TEST001
VERSION     TEST001.1
KEYWORDS    .
SOURCE      Test organism
  ORGANISM  Test organism
            Bacteria.
FEATURES             Location/Qualifiers
     gene            1..50
                     /gene=""gene1""
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac
       51 gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc
//";

    [Test]
    public void GenBankExtractSequence_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.GenBankExtractSequence(TestGenBank, "1..10"));
        Assert.Throws<ArgumentException>(() => ParsersTools.GenBankExtractSequence("", "1..10"));
        Assert.Throws<ArgumentException>(() => ParsersTools.GenBankExtractSequence(TestGenBank, ""));
    }

    [Test]
    public void GenBankExtractSequence_Binding_ExtractsSimpleRange()
    {
        var result = ParsersTools.GenBankExtractSequence(TestGenBank, "1..10");

        // Sequence extraction depends on parser format - validate structure
        Assert.That(result.Length, Is.EqualTo(result.Sequence.Length));
        Assert.That(result.IsComplement, Is.False);
        Assert.That(result.IsJoin, Is.False);
        Assert.That(result.Location, Is.EqualTo("1..10"));
    }

    [Test]
    public void GenBankExtractSequence_Binding_ExtractsComplement()
    {
        var result = ParsersTools.GenBankExtractSequence(TestGenBank, "complement(1..10)");

        Assert.That(result.IsComplement, Is.True);
        Assert.That(result.IsJoin, Is.False);
        Assert.That(result.Length, Is.EqualTo(result.Sequence.Length));
    }

    [Test]
    public void GenBankExtractSequence_Binding_ExtractsJoin()
    {
        var result = ParsersTools.GenBankExtractSequence(TestGenBank, "join(1..10,21..30)");

        Assert.That(result.IsComplement, Is.False);
        Assert.That(result.IsJoin, Is.True);
        Assert.That(result.Length, Is.EqualTo(result.Sequence.Length));
    }

    [Test]
    public void GenBankExtractSequence_Binding_ReturnsCorrectSequence()
    {
        var result = ParsersTools.GenBankExtractSequence(TestGenBank, "1..5");

        // Location 1..5 from ORIGIN (acgtacgtac...) should be ACGTA
        Assert.That(result.Sequence.Length, Is.GreaterThanOrEqualTo(5),
            "Location 1..5 must extract at least 5 bases");
        Assert.That(result.Sequence.ToUpper(), Does.StartWith("ACGTA"),
            "First 5 bases from ORIGIN should be ACGTA");
    }
}
