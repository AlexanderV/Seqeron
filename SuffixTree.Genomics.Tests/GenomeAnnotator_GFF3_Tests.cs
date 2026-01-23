using NUnit.Framework;
using System.Collections.Generic;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

/// <summary>
/// Evidence-based tests for GenomeAnnotator GFF3 I/O methods.
/// Test Unit: ANNOT-GFF-001
/// 
/// Sources:
/// - Sequence Ontology GFF3 Specification v1.26 (Lincoln Stein, 2020)
/// - Wikipedia: General feature format
/// - RFC 3986: URI percent-encoding
/// </summary>
[TestFixture]
public class GenomeAnnotator_GFF3_Tests
{
    #region ParseGff3 - Basic Parsing

    /// <summary>
    /// M1: ParseGff3 parses valid GFF3 line correctly.
    /// Source: GFF3 Specification - 9 tab-delimited columns
    /// </summary>
    [Test]
    public void ParseGff3_ValidLine_ParsesCorrectly()
    {
        var lines = new List<string>
        {
            "##gff-version 3",
            "seq1\t.\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=testGene"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(features[0].Type, Is.EqualTo("gene"));
            Assert.That(features[0].Start, Is.EqualTo(100));
            Assert.That(features[0].End, Is.EqualTo(500));
            Assert.That(features[0].Strand, Is.EqualTo('+'));
        });
    }

    /// <summary>
    /// M2: ParseGff3 extracts all 9 columns correctly.
    /// Source: GFF3 Specification - Column definitions
    /// </summary>
    [Test]
    public void ParseGff3_AllColumns_ExtractsCorrectly()
    {
        var lines = new List<string>
        {
            "chr1\tENSEMBL\tCDS\t1000\t2000\t95.5\t-\t2\tID=cds1;Name=TestCDS"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features, Has.Count.EqualTo(1));
        var f = features[0];
        Assert.Multiple(() =>
        {
            Assert.That(f.FeatureId, Is.EqualTo("cds1"), "Column 9: ID attribute");
            Assert.That(f.Type, Is.EqualTo("CDS"), "Column 3: type");
            Assert.That(f.Start, Is.EqualTo(1000), "Column 4: start");
            Assert.That(f.End, Is.EqualTo(2000), "Column 5: end");
            Assert.That(f.Score, Is.EqualTo(95.5).Within(0.01), "Column 6: score");
            Assert.That(f.Strand, Is.EqualTo('-'), "Column 7: strand");
            Assert.That(f.Phase, Is.EqualTo(2), "Column 8: phase");
        });
    }

    /// <summary>
    /// M18: ParseGff3 handles all strand values (+, -, .).
    /// Source: GFF3 Specification - Column 7 strand values
    /// </summary>
    [TestCase('+', '+')]
    [TestCase('-', '-')]
    [TestCase('.', '.')]
    public void ParseGff3_StrandValues_ParsedCorrectly(char inputStrand, char expectedStrand)
    {
        var lines = new List<string>
        {
            $"seq1\t.\tgene\t100\t500\t.\t{inputStrand}\t.\tID=g1"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features[0].Strand, Is.EqualTo(expectedStrand));
    }

    #endregion

    #region ParseGff3 - Score Handling

    /// <summary>
    /// M3: ParseGff3 parses numeric score correctly.
    /// Source: GFF3 Specification - Column 6 score as floating point
    /// </summary>
    [Test]
    public void ParseGff3_WithScore_ParsesScore()
    {
        var lines = new List<string>
        {
            "seq1\t.\tCDS\t100\t500\t0.95\t+\t0\tID=cds1"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features[0].Score, Is.EqualTo(0.95).Within(0.001));
    }

    /// <summary>
    /// M4: ParseGff3 handles "." as null score.
    /// Source: GFF3 Specification - "." denotes undefined value
    /// </summary>
    [Test]
    public void ParseGff3_UndefinedScore_ReturnsNull()
    {
        var lines = new List<string>
        {
            "seq1\t.\tgene\t100\t500\t.\t+\t.\tID=g1"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features[0].Score, Is.Null);
    }

    #endregion

    #region ParseGff3 - Phase Handling

    /// <summary>
    /// M5: ParseGff3 parses phase as integer (0, 1, 2).
    /// Source: GFF3 Specification - Column 8 phase for CDS features
    /// </summary>
    [TestCase("0", 0)]
    [TestCase("1", 1)]
    [TestCase("2", 2)]
    public void ParseGff3_PhaseValues_ParsedAsInteger(string phaseInput, int expectedPhase)
    {
        var lines = new List<string>
        {
            $"seq1\t.\tCDS\t100\t500\t.\t+\t{phaseInput}\tID=cds1"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features[0].Phase, Is.EqualTo(expectedPhase));
    }

    /// <summary>
    /// M6: ParseGff3 handles "." as null phase.
    /// Source: GFF3 Specification - "." for non-CDS features
    /// </summary>
    [Test]
    public void ParseGff3_UndefinedPhase_ReturnsNull()
    {
        var lines = new List<string>
        {
            "seq1\t.\tgene\t100\t500\t.\t+\t.\tID=g1"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features[0].Phase, Is.Null);
    }

    #endregion

    #region ParseGff3 - Comment and Directive Handling

    /// <summary>
    /// M7: ParseGff3 skips comment lines starting with single #.
    /// Source: GFF3 Specification - Lines beginning with # are comments
    /// </summary>
    [Test]
    public void ParseGff3_SkipsComments()
    {
        var lines = new List<string>
        {
            "# This is a comment",
            "##gff-version 3",
            "seq1\t.\tgene\t100\t500\t.\t+\t.\tID=gene1"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features, Has.Count.EqualTo(1));
    }

    /// <summary>
    /// M8: ParseGff3 skips directive lines starting with ##.
    /// Source: GFF3 Specification - ## denotes meta-directives
    /// </summary>
    [Test]
    public void ParseGff3_SkipsDirectives()
    {
        var lines = new List<string>
        {
            "##gff-version 3",
            "##sequence-region seq1 1 10000",
            "##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606",
            "seq1\t.\tgene\t100\t500\t.\t+\t.\tID=gene1"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features, Has.Count.EqualTo(1));
        Assert.That(features[0].FeatureId, Is.EqualTo("gene1"));
    }

    /// <summary>
    /// M9: ParseGff3 skips empty lines.
    /// Source: GFF3 Specification - Empty lines should be ignored
    /// </summary>
    [Test]
    public void ParseGff3_SkipsEmptyLines()
    {
        var lines = new List<string>
        {
            "",
            "   ",
            "seq1\t.\tgene\t100\t500\t.\t+\t.\tID=gene1",
            "",
            "seq1\t.\tgene\t600\t900\t.\t+\t.\tID=gene2"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features, Has.Count.EqualTo(2));
    }

    #endregion

    #region ParseGff3 - Attribute Parsing

    /// <summary>
    /// M10: ParseGff3 parses semicolon-separated attributes.
    /// Source: GFF3 Specification - Column 9 attributes as key=value pairs
    /// </summary>
    [Test]
    public void ParseGff3_ParsesAttributes()
    {
        var lines = new List<string>
        {
            "seq1\t.\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=myGene;product=test%20protein"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(features[0].Attributes["ID"], Is.EqualTo("gene1"));
            Assert.That(features[0].Attributes["Name"], Is.EqualTo("myGene"));
            Assert.That(features[0].Attributes["product"], Is.EqualTo("test protein"));
        });
    }

    /// <summary>
    /// M11: ParseGff3 decodes URL-encoded attribute values.
    /// Source: RFC 3986 - Percent-encoding
    /// </summary>
    [TestCase("product=test%20protein", "product", "test protein")]
    [TestCase("Name=gene%3B1", "Name", "gene;1")]
    [TestCase("Note=equals%3Dtest", "Note", "equals=test")]
    public void ParseGff3_DecodesUrlEncodedValues(string attrString, string key, string expectedValue)
    {
        var lines = new List<string>
        {
            $"seq1\t.\tgene\t100\t500\t.\t+\t.\tID=g1;{attrString}"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features[0].Attributes[key], Is.EqualTo(expectedValue));
    }

    /// <summary>
    /// S1: ParseGff3 handles multiple attributes on single line.
    /// </summary>
    [Test]
    public void ParseGff3_MultipleAttributes_AllParsed()
    {
        var lines = new List<string>
        {
            "seq1\t.\tCDS\t100\t500\t.\t+\t0\tID=cds1;Name=CDS1;Parent=mRNA1;product=hypothetical;Note=test"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features[0].Attributes, Has.Count.GreaterThanOrEqualTo(4));
    }

    #endregion

    #region ParseGff3 - Error Handling

    /// <summary>
    /// M12: ParseGff3 skips malformed lines with fewer than 9 fields.
    /// Source: GFF3 Specification - 9 columns required
    /// </summary>
    [Test]
    public void ParseGff3_MalformedLine_Skipped()
    {
        var lines = new List<string>
        {
            "seq1\t.\tgene",  // Only 3 fields
            "seq1\t.\tgene\t100\t500",  // Only 5 fields
            "seq1\t.\tgene\t100\t500\t.\t+\t.\tID=valid"  // 9 fields - valid
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features, Has.Count.EqualTo(1));
        Assert.That(features[0].FeatureId, Is.EqualTo("valid"));
    }

    /// <summary>
    /// M13: ParseGff3 auto-generates ID if ID attribute is missing.
    /// Source: Implementation behavior
    /// </summary>
    [Test]
    public void ParseGff3_MissingIdAttribute_AutoGeneratesId()
    {
        var lines = new List<string>
        {
            "seq1\t.\tgene\t100\t500\t.\t+\t.\tName=testGene"  // No ID attribute
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features[0].FeatureId, Is.Not.Null.And.Not.Empty);
        Assert.That(features[0].FeatureId, Does.StartWith("feature_"));
    }

    /// <summary>
    /// Empty input returns empty result.
    /// Source: GFF3 Specification
    /// </summary>
    [Test]
    public void ParseGff3_EmptyInput_ReturnsEmpty()
    {
        var lines = new List<string>();

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features, Is.Empty);
    }

    /// <summary>
    /// Input with only comments returns empty result.
    /// </summary>
    [Test]
    public void ParseGff3_OnlyComments_ReturnsEmpty()
    {
        var lines = new List<string>
        {
            "##gff-version 3",
            "# Comment 1",
            "# Comment 2"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features, Is.Empty);
    }

    #endregion

    #region ToGff3 - Basic Export

    /// <summary>
    /// M14: ToGff3 emits ##gff-version 3 header.
    /// Source: GFF3 Specification - Version directive required
    /// </summary>
    [Test]
    public void ToGff3_EmitsVersionHeader()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new(
                GeneId: "gene1",
                Start: 99,
                End: 500,
                Strand: '+',
                Type: "CDS",
                Product: "test",
                Attributes: new Dictionary<string, string>())
        };

        var lines = GenomeAnnotator.ToGff3(annotations, "chr1").ToList();

        Assert.That(lines[0], Is.EqualTo("##gff-version 3"));
    }

    /// <summary>
    /// M15: ToGff3 converts 0-based internal Start to 1-based GFF3 coordinate.
    /// Source: GFF3 Specification - 1-based coordinates
    /// </summary>
    [Test]
    public void ToGff3_ConvertsToOneBased()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new(
                GeneId: "gene1",
                Start: 99,   // 0-based internal
                End: 500,
                Strand: '+',
                Type: "CDS",
                Product: "test",
                Attributes: new Dictionary<string, string>())
        };

        var lines = GenomeAnnotator.ToGff3(annotations, "chr1").ToList();
        var dataLine = lines[1];
        var fields = dataLine.Split('\t');

        Assert.That(fields[3], Is.EqualTo("100"), "Start should be 1-based (99 + 1 = 100)");
    }

    /// <summary>
    /// M17: ToGff3 produces valid tab-delimited output with 9 columns.
    /// Source: GFF3 Specification - 9 tab-delimited columns
    /// </summary>
    [Test]
    public void ToGff3_GeneratesValidTabDelimitedOutput()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new(
                GeneId: "gene1",
                Start: 99,
                End: 500,
                Strand: '+',
                Type: "CDS",
                Product: "hypothetical protein",
                Attributes: new Dictionary<string, string> { ["frame"] = "1" })
        };

        var lines = GenomeAnnotator.ToGff3(annotations, "chr1").ToList();

        Assert.That(lines, Has.Count.EqualTo(2));  // Header + 1 feature
        Assert.That(lines[0], Does.Contain("##gff-version 3"));

        var dataLine = lines[1];
        var fields = dataLine.Split('\t');
        Assert.That(fields, Has.Length.EqualTo(9), "Should have exactly 9 tab-separated fields");
        Assert.Multiple(() =>
        {
            Assert.That(fields[0], Is.EqualTo("chr1"), "Column 1: seqid");
            Assert.That(fields[2], Is.EqualTo("CDS"), "Column 3: type");
            Assert.That(fields[8], Does.Contain("ID=gene1"), "Column 9: attributes");
        });
    }

    #endregion

    #region ToGff3 - Encoding

    /// <summary>
    /// M16: ToGff3 URL-encodes special characters in attribute values.
    /// Source: RFC 3986 - Percent-encoding for special characters
    /// </summary>
    [Test]
    public void ToGff3_EscapesSpecialCharacters()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new(
                GeneId: "gene 1",  // Space in ID
                Start: 0,
                End: 100,
                Strand: '+',
                Type: "CDS",
                Product: "test;product",  // Semicolon in product
                Attributes: new Dictionary<string, string>())
        };

        var lines = GenomeAnnotator.ToGff3(annotations).ToList();
        var dataLine = lines[1];

        Assert.That(dataLine, Does.Contain("gene%201"), "Space should be encoded as %20");
        Assert.That(dataLine, Does.Contain("test%3Bproduct"), "Semicolon should be encoded as %3B");
    }

    /// <summary>
    /// S3: ToGff3 excludes translation attribute to avoid bloat.
    /// Source: Implementation design - translation can be very large
    /// </summary>
    [Test]
    public void ToGff3_ExcludesTranslationAttribute()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new(
                GeneId: "gene1",
                Start: 0,
                End: 100,
                Strand: '+',
                Type: "CDS",
                Product: "test",
                Attributes: new Dictionary<string, string>
                {
                    ["translation"] = "MKKKKKKKKKKKKKKKKKKKKKKKKKKK",
                    ["Note"] = "important"
                })
        };

        var lines = GenomeAnnotator.ToGff3(annotations).ToList();
        var dataLine = lines[1];

        Assert.Multiple(() =>
        {
            Assert.That(dataLine, Does.Not.Contain("translation="));
            Assert.That(dataLine, Does.Contain("Note=important"));
        });
    }

    #endregion

    #region ToGff3 - Edge Cases

    /// <summary>
    /// C2: ToGff3 handles empty annotations list.
    /// </summary>
    [Test]
    public void ToGff3_EmptyAnnotations_ReturnsOnlyHeader()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>();

        var lines = GenomeAnnotator.ToGff3(annotations).ToList();

        Assert.That(lines, Has.Count.EqualTo(1));
        Assert.That(lines[0], Is.EqualTo("##gff-version 3"));
    }

    /// <summary>
    /// ToGff3 uses default seqId when not specified.
    /// </summary>
    [Test]
    public void ToGff3_DefaultSeqId_UsesSeq1()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new(
                GeneId: "gene1",
                Start: 0,
                End: 100,
                Strand: '+',
                Type: "gene",
                Product: "test",
                Attributes: new Dictionary<string, string>())
        };

        var lines = GenomeAnnotator.ToGff3(annotations).ToList();
        var fields = lines[1].Split('\t');

        Assert.That(fields[0], Is.EqualTo("seq1"));
    }

    #endregion

    #region Roundtrip Tests

    /// <summary>
    /// S4: Parse(ToGff3()) preserves core feature data.
    /// Note: Not exact roundtrip due to different record types, but core data preserved.
    /// </summary>
    [Test]
    public void Roundtrip_CoreDataPreserved()
    {
        // Create annotation and export
        var original = new GenomeAnnotator.GeneAnnotation(
            GeneId: "gene1",
            Start: 99,  // 0-based
            End: 500,
            Strand: '+',
            Type: "CDS",
            Product: "test protein",
            Attributes: new Dictionary<string, string> { ["Note"] = "important" });

        var exported = GenomeAnnotator.ToGff3(new[] { original }, "chr1").ToList();

        // Parse back
        var parsed = GenomeAnnotator.ParseGff3(exported).ToList();

        Assert.That(parsed, Has.Count.EqualTo(1));
        var f = parsed[0];
        Assert.Multiple(() =>
        {
            Assert.That(f.FeatureId, Is.EqualTo("gene1"));
            Assert.That(f.Start, Is.EqualTo(100), "1-based after roundtrip");
            Assert.That(f.End, Is.EqualTo(500));
            Assert.That(f.Strand, Is.EqualTo('+'));
            Assert.That(f.Type, Is.EqualTo("CDS"));
        });
    }

    #endregion

    #region Multiple Features

    /// <summary>
    /// ParseGff3 handles multiple features in correct order.
    /// </summary>
    [Test]
    public void ParseGff3_MultipleFeatures_PreservesOrder()
    {
        var lines = new List<string>
        {
            "##gff-version 3",
            "chr1\t.\tgene\t100\t500\t.\t+\t.\tID=gene1",
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=mrna1;Parent=gene1",
            "chr1\t.\texon\t100\t200\t.\t+\t.\tID=exon1;Parent=mrna1",
            "chr1\t.\texon\t300\t500\t.\t+\t.\tID=exon2;Parent=mrna1"
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();

        Assert.That(features, Has.Count.EqualTo(4));
        Assert.Multiple(() =>
        {
            Assert.That(features[0].Type, Is.EqualTo("gene"));
            Assert.That(features[1].Type, Is.EqualTo("mRNA"));
            Assert.That(features[2].Type, Is.EqualTo("exon"));
            Assert.That(features[3].Type, Is.EqualTo("exon"));
        });
    }

    #endregion
}
