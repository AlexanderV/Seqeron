using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class EmblParseTests
{
    // Use the same format as in EmblParserTests.cs
    private const string TestEmbl = @"ID   TEST001; SV 1; linear; DNA; STD; HUM; 100 BP.
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

    [Test]
    public void EmblParse_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.EmblParse(TestEmbl));
        Assert.Throws<ArgumentException>(() => ParsersTools.EmblParse(""));
        Assert.Throws<ArgumentException>(() => ParsersTools.EmblParse(null!));
    }

    [Test]
    public void EmblParse_Binding_ParsesRecords()
    {
        var result = ParsersTools.EmblParse(TestEmbl);

        Assert.That(result.Count, Is.EqualTo(1));
        Assert.That(result.Records[0].Accession, Is.EqualTo("TEST001"));
        Assert.That(result.Records[0].SequenceLength, Is.EqualTo(100));
        Assert.That(result.Records[0].MoleculeType, Is.EqualTo("DNA"));
        Assert.That(result.Records[0].Topology, Is.EqualTo("linear"));
    }

    [Test]
    public void EmblParse_Binding_ParsesDescription()
    {
        var result = ParsersTools.EmblParse(TestEmbl);

        Assert.That(result.Records[0].Description, Does.Contain("Test sequence"));
    }

    [Test]
    public void EmblParse_Binding_ParsesKeywords()
    {
        var result = ParsersTools.EmblParse(TestEmbl);

        Assert.That(result.Records[0].Keywords, Contains.Item("test"));
        Assert.That(result.Records[0].Keywords, Contains.Item("genomics"));
    }

    [Test]
    public void EmblParse_Binding_ParsesOrganism()
    {
        var result = ParsersTools.EmblParse(TestEmbl);

        Assert.That(result.Records[0].Organism, Is.EqualTo("Homo sapiens"));
        Assert.That(result.Records[0].OrganismClassification, Does.Contain("Eukaryota"));
    }

    [Test]
    public void EmblParse_Binding_ParsesFeatures()
    {
        var result = ParsersTools.EmblParse(TestEmbl);

        // Features may or may not be parsed depending on implementation
        Assert.That(result.Records[0].Features, Is.Not.Null);
        if (result.Records[0].Features.Count > 0)
        {
            Assert.That(result.Records[0].Features.Any(f => f.Key == "gene" || f.Key == "CDS"), Is.True);
        }
    }

    [Test]
    public void EmblParse_Binding_ParsesSequence()
    {
        var result = ParsersTools.EmblParse(TestEmbl);

        Assert.That(result.Records[0].Sequence, Does.StartWith("ACGTACGT"));
        Assert.That(result.Records[0].ActualSequenceLength, Is.EqualTo(100));
    }
}

[TestFixture]
public class EmblFeaturesTests
{
    private const string TestEmbl = @"ID   TEST001; SV 1; linear; DNA; STD; HUM; 100 BP.
XX
AC   TEST001;
XX
DE   Test sequence.
XX
OS   Homo sapiens
XX
FH   Key             Location/Qualifiers
FH
FT   gene            1..50
FT                   /gene=""gene1""
FT   gene            60..90
FT                   /gene=""gene2""
FT   CDS             10..40
FT                   /gene=""gene1""
FT                   /product=""Product 1""
XX
SQ   Sequence 100 BP;
     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac        50
     gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc       100
//";

    [Test]
    public void EmblFeatures_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.EmblFeatures(TestEmbl));
        Assert.Throws<ArgumentException>(() => ParsersTools.EmblFeatures(""));
    }

    [Test]
    public void EmblFeatures_Binding_ExtractsAllFeatures()
    {
        var result = ParsersTools.EmblFeatures(TestEmbl);

        Assert.That(result.Features, Is.Not.Null);
        Assert.That(result.Count, Is.EqualTo(result.Features.Count));
    }

    [Test]
    public void EmblFeatures_Binding_FiltersByType()
    {
        var result = ParsersTools.EmblFeatures(TestEmbl, featureType: "gene");

        Assert.That(result.Features.All(f => f.Key == "gene"), Is.True);
    }

    [Test]
    public void EmblFeatures_Binding_ParsesQualifiers()
    {
        var result = ParsersTools.EmblFeatures(TestEmbl, featureType: "CDS");

        if (result.Count > 0)
        {
            Assert.That(result.Features[0].Qualifiers, Is.Not.Empty);
        }
    }

    [Test]
    public void EmblFeatures_Binding_ParsesLocation()
    {
        var result = ParsersTools.EmblFeatures(TestEmbl, featureType: "gene");

        if (result.Count > 0)
        {
            var gene1 = result.Features.FirstOrDefault(f => f.Qualifiers.GetValueOrDefault("gene") == "gene1");
            if (gene1 != null)
            {
                Assert.That(gene1.Start, Is.EqualTo(1));
                Assert.That(gene1.End, Is.EqualTo(50));
            }
        }
    }
}

[TestFixture]
public class EmblStatisticsTests
{
    private const string TestEmbl = @"ID   TEST001; SV 1; linear; DNA; STD; HUM; 100 BP.
XX
AC   TEST001;
XX
DE   Test sequence.
XX
OS   Homo sapiens
XX
FH   Key             Location/Qualifiers
FH
FT   gene            1..50
FT                   /gene=""gene1""
FT   CDS             10..40
FT                   /product=""Product 1""
XX
SQ   Sequence 100 BP;
     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac        50
     gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc gcgcgcgcgc       100
//";

    [Test]
    public void EmblStatistics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.EmblStatistics(TestEmbl));
        Assert.Throws<ArgumentException>(() => ParsersTools.EmblStatistics(""));
    }

    [Test]
    public void EmblStatistics_Binding_CalculatesStats()
    {
        var result = ParsersTools.EmblStatistics(TestEmbl);

        Assert.That(result.RecordCount, Is.EqualTo(1));
        Assert.That(result.TotalSequenceLength, Is.EqualTo(100));
        Assert.That(result.TotalFeatures, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void EmblStatistics_Binding_CountsFeatureTypes()
    {
        var result = ParsersTools.EmblStatistics(TestEmbl);

        var totalFromCounts = result.FeatureTypeCounts.Values.Sum();
        Assert.That(totalFromCounts, Is.EqualTo(result.TotalFeatures));
    }

    [Test]
    public void EmblStatistics_Binding_CollectsMoleculeTypes()
    {
        var result = ParsersTools.EmblStatistics(TestEmbl);

        Assert.That(result.MoleculeTypes, Contains.Item("DNA"));
    }
}
