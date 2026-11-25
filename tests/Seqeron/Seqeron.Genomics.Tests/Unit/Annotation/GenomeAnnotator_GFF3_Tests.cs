namespace Seqeron.Genomics.Tests.Unit.Annotation;

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
    /// M18: ParseGff3 handles all strand values (+, -, ., ?).
    /// Source: GFF3 Specification - Column 7 strand values
    /// "? can be used for features whose strandedness is relevant, but unknown."
    /// </summary>
    [TestCase('+', '+')]
    [TestCase('-', '-')]
    [TestCase('.', '.')]
    [TestCase('?', '?')]
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
        Assert.That(features[0].FeatureId, Is.EqualTo("gene1"));
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
        Assert.Multiple(() =>
        {
            Assert.That(features[0].FeatureId, Is.EqualTo("gene1"));
            Assert.That(features[1].FeatureId, Is.EqualTo("gene2"));
        });
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
    [TestCase("Name=a%09b", "Name", "a\tb")]
    [TestCase("Note=line1%0Aline2", "Note", "line1\nline2")]
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

        Assert.That(features[0].Attributes, Has.Count.EqualTo(5));
        Assert.Multiple(() =>
        {
            Assert.That(features[0].Attributes["ID"], Is.EqualTo("cds1"));
            Assert.That(features[0].Attributes["Name"], Is.EqualTo("CDS1"));
            Assert.That(features[0].Attributes["Parent"], Is.EqualTo("mRNA1"));
            Assert.That(features[0].Attributes["product"], Is.EqualTo("hypothetical"));
            Assert.That(features[0].Attributes["Note"], Is.EqualTo("test"));
        });
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

        var dataLine = lines[1];
        var fields = dataLine.Split('\t');
        Assert.That(fields, Has.Length.EqualTo(9), "Should have exactly 9 tab-separated fields");
        Assert.Multiple(() =>
        {
            Assert.That(fields[0], Is.EqualTo("chr1"), "Column 1: seqid");
            Assert.That(fields[1], Is.EqualTo("."), "Column 2: source");
            Assert.That(fields[2], Is.EqualTo("CDS"), "Column 3: type");
            Assert.That(fields[3], Is.EqualTo("100"), "Column 4: start (0-based 99 + 1)");
            Assert.That(fields[4], Is.EqualTo("500"), "Column 5: end");
            Assert.That(fields[5], Is.EqualTo("."), "Column 6: score");
            Assert.That(fields[6], Is.EqualTo("+"), "Column 7: strand");
            Assert.That(fields[7], Is.EqualTo("0"), "Column 8: phase (CDS)");
            Assert.That(fields[8], Does.Contain("ID=gene1"), "Column 9: attributes contain ID");
            Assert.That(fields[8], Does.Contain("product=hypothetical protein"), "Column 9: attributes contain product");
        });
    }

    #endregion

    #region ToGff3 - Encoding

    /// <summary>
    /// M16: ToGff3 encodes only GFF3-required special characters in attribute values.
    /// Source: GFF3 Spec v1.26 — semicolons, equals, ampersands, commas MUST be encoded;
    /// "no other characters may be encoded"; "unescaped spaces are allowed within fields."
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

        Assert.Multiple(() =>
        {
            Assert.That(dataLine, Does.Contain("ID=gene 1"),
                "Spaces must not be encoded per GFF3 spec: 'unescaped spaces are allowed within fields'");
            Assert.That(dataLine, Does.Not.Contain("gene%201"),
                "Spaces must not be percent-encoded per GFF3 spec");
            Assert.That(dataLine, Does.Contain("test%3Bproduct"),
                "Semicolons must be encoded as %3B per GFF3 spec");
        });
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

    /// <summary>
    /// M19: ToGff3 encodes all GFF3-required column 9 characters correctly.
    /// Source: GFF3 Spec v1.26 — in column 9, semicolon (%3B), equals (%3D),
    /// ampersand (%26), comma (%2C) must be encoded.
    /// </summary>
    [Test]
    public void ToGff3_EncodesAllRequiredGff3Characters()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new(
                GeneId: "g1",
                Start: 0,
                End: 100,
                Strand: '+',
                Type: "gene",
                Product: "a;b=c&d,e",
                Attributes: new Dictionary<string, string>())
        };

        var lines = GenomeAnnotator.ToGff3(annotations).ToList();
        var dataLine = lines[1];

        Assert.Multiple(() =>
        {
            Assert.That(dataLine, Does.Contain("a%3Bb%3Dc%26d%2Ce"),
                "Semicolon, equals, ampersand, comma must all be encoded");
        });
    }

    /// <summary>
    /// M23: ToGff3 emits the feature's real source (column 2) when present.
    /// Source: SO GFF3 Spec v1.26 — column 2 is "a free text qualifier intended to describe
    /// the algorithm or operating procedure that generated this feature"; the "." placeholder
    /// is only for undefined fields.
    /// </summary>
    [Test]
    public void ToGff3_WithSource_EmitsRealSource()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new(
                GeneId: "gene1",
                Start: 99,
                End: 500,
                Strand: '+',
                Type: "gene",
                Product: "test",
                Attributes: new Dictionary<string, string>(),
                Source: "ENSEMBL")
        };

        var fields = GenomeAnnotator.ToGff3(annotations, "chr1").ToList()[1].Split('\t');

        Assert.That(fields[1], Is.EqualTo("ENSEMBL"),
            "Column 2 must carry the feature's real source when present (GFF3 Spec v1.26 column 2)");
    }

    /// <summary>
    /// M24: ToGff3 falls back to "." for source when none was supplied.
    /// Source: SO GFF3 Spec v1.26 — "Undefined fields are replaced with the '.' character."
    /// </summary>
    [Test]
    public void ToGff3_NoSource_EmitsDotPlaceholder()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new(
                GeneId: "gene1",
                Start: 99,
                End: 500,
                Strand: '+',
                Type: "gene",
                Product: "test",
                Attributes: new Dictionary<string, string>())
        };

        var fields = GenomeAnnotator.ToGff3(annotations, "chr1").ToList()[1].Split('\t');

        Assert.That(fields[1], Is.EqualTo("."),
            "Column 2 must be the '.' placeholder when source is absent (GFF3 Spec v1.26)");
    }

    /// <summary>
    /// M25: ToGff3 emits the feature's real score (column 6) when present, culture-invariant.
    /// Source: SO GFF3 Spec v1.26 — column 6 is "the score of the feature, a floating point number".
    /// </summary>
    [Test]
    public void ToGff3_WithScore_EmitsRealScore()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new(
                GeneId: "gene1",
                Start: 99,
                End: 500,
                Strand: '+',
                Type: "gene",
                Product: "test",
                Attributes: new Dictionary<string, string>(),
                Source: ".",
                Score: 95.5)
        };

        var fields = GenomeAnnotator.ToGff3(annotations, "chr1").ToList()[1].Split('\t');

        Assert.That(fields[5], Is.EqualTo("95.5"),
            "Column 6 must carry the feature's real floating-point score (GFF3 Spec v1.26 column 6)");
    }

    /// <summary>
    /// M26: ToGff3 falls back to "." for score when none was supplied.
    /// Source: SO GFF3 Spec v1.26 — "Undefined fields are replaced with the '.' character."
    /// </summary>
    [Test]
    public void ToGff3_NoScore_EmitsDotPlaceholder()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new(
                GeneId: "gene1",
                Start: 99,
                End: 500,
                Strand: '+',
                Type: "gene",
                Product: "test",
                Attributes: new Dictionary<string, string>())
        };

        var fields = GenomeAnnotator.ToGff3(annotations, "chr1").ToList()[1].Split('\t');

        Assert.That(fields[5], Is.EqualTo("."),
            "Column 6 must be the '.' placeholder when score is absent (GFF3 Spec v1.26)");
    }

    /// <summary>
    /// M27: ToGff3 round-trips real source and score through ParseGff3.
    /// Source: SO GFF3 Spec v1.26 — columns 2 and 6 carry real values.
    /// </summary>
    [Test]
    public void ToGff3_SourceAndScore_SurviveRoundtrip()
    {
        var original = new GenomeAnnotator.GeneAnnotation(
            GeneId: "gene1",
            Start: 99,
            End: 500,
            Strand: '+',
            Type: "gene",
            Product: "test protein",
            Attributes: new Dictionary<string, string>(),
            Source: "ENSEMBL",
            Score: 12.75);

        var exported = GenomeAnnotator.ToGff3(new[] { original }, "chr1").ToList();
        var parsed = GenomeAnnotator.ParseGff3(exported).ToList();

        Assert.That(parsed, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            // Parser exposes source as parts[1] only via raw line; verify via the emitted line and Score.
            Assert.That(exported[1].Split('\t')[1], Is.EqualTo("ENSEMBL"), "source survives export");
            Assert.That(parsed[0].Score, Is.EqualTo(12.75).Within(1e-10), "score survives the round-trip");
        });
    }

    /// <summary>
    /// M28: ToGff3 computes cumulative CDS phase across a multi-segment transcript on the
    /// PLUS strand. Hand-derived from the SO GFF3 Spec v1.26 canonical-gene example
    /// (mRNA00003 / cds00003): three + strand CDS segments of lengths 602, 501, 602
    /// (1-based 3301-3902, 5000-5500, 7000-7600) yield phases 0, 1, 1.
    ///   seg1: (3 − (0   mod 3)) mod 3 = 0
    ///   seg2: (3 − (602 mod 3)) mod 3 = (3 − 2) mod 3 = 1
    ///   seg3: (3 − (1103 mod 3)) mod 3 = (3 − 2) mod 3 = 1
    /// Source: SO GFF3 Spec v1.26 — column 8 phase definition + canonical-gene example.
    /// </summary>
    [Test]
    public void ToGff3_PlusStrandMultiSegmentCds_ComputesCumulativePhase()
    {
        // Internal 0-based half-open: 1-based start S => internal Start = S-1; End = 1-based end.
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new("mRNA3", 3300, 3902, '+', "CDS", "edenprotein", new Dictionary<string, string>()), // 1-based 3301-3902, len 602
            new("mRNA3", 4999, 5500, '+', "CDS", "edenprotein", new Dictionary<string, string>()), // 1-based 5000-5500, len 501
            new("mRNA3", 6999, 7600, '+', "CDS", "edenprotein", new Dictionary<string, string>()), // 1-based 7000-7600, len 602
        };

        var lines = GenomeAnnotator.ToGff3(annotations, "ctg123").ToList();
        var phases = lines.Skip(1).Select(l => l.Split('\t')[7]).ToList();

        Assert.That(phases, Is.EqualTo(new[] { "0", "1", "1" }),
            "Plus-strand CDS phases must match the SO canonical-gene example (0, 1, 1)");
    }

    /// <summary>
    /// M29: ToGff3 computes cumulative CDS phase across a multi-segment transcript on the
    /// MINUS strand, where the 5' end is the feature's END, so segments are ordered by
    /// DESCENDING genomic coordinate. Hand-derived from SO GFF3 Spec v1.26 column-8 rules:
    /// three − strand CDS segments (1-based 7000-7600 len 601, 5000-5500 len 501, 3301-3902 len 602)
    /// processed 5'→3' (highest coord first):
    ///   seg @7000 (5'-most): (3 − (0   mod 3)) mod 3 = 0
    ///   seg @5000:           (3 − (601 mod 3)) mod 3 = (3 − 1) mod 3 = 2   [601 mod 3 = 1]
    ///   seg @3301:           (3 − (1102 mod 3)) mod 3 = (3 − 1) mod 3 = 2   [1102 mod 3 = 1]
    /// Emitted in input (ascending genomic) order 3301, 5000, 7000 → phases 2, 2, 0.
    /// Source: SO GFF3 Spec v1.26 — "the 5' end for CDS features on the minus strand is the
    /// feature's end"; phase = bases forward to the next codon.
    /// </summary>
    [Test]
    public void ToGff3_MinusStrandMultiSegmentCds_ComputesCumulativePhase()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new("mRNAneg", 3300, 3902, '-', "CDS", "p", new Dictionary<string, string>()), // 1-based 3301-3902, len 602
            new("mRNAneg", 4999, 5500, '-', "CDS", "p", new Dictionary<string, string>()), // 1-based 5000-5500, len 501
            new("mRNAneg", 6999, 7600, '-', "CDS", "p", new Dictionary<string, string>()), // 1-based 7000-7600, len 601
        };

        var lines = GenomeAnnotator.ToGff3(annotations, "ctg123").ToList();
        var phases = lines.Skip(1).Select(l => l.Split('\t')[7]).ToList();

        Assert.That(phases, Is.EqualTo(new[] { "2", "2", "0" }),
            "Minus-strand CDS phases (input order 3301,5000,7000) must be 2, 2, 0 with the 5'-most (7000) segment at phase 0");
    }

    /// <summary>
    /// M30: ToGff3 keeps separate transcripts independent — phase accumulation never bleeds
    /// across different GeneIds. Two single-segment CDS transcripts each begin at phase 0.
    /// Source: SO GFF3 Spec v1.26 — phase is relative to the start of the current CDS feature
    /// within its own transcript.
    /// </summary>
    [Test]
    public void ToGff3_DistinctTranscripts_PhaseDoesNotAccumulateAcross()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new("txA", 0, 100, '+', "CDS", "a", new Dictionary<string, string>()),   // len 100
            new("txB", 200, 305, '+', "CDS", "b", new Dictionary<string, string>()), // separate transcript
        };

        var lines = GenomeAnnotator.ToGff3(annotations, "chr1").ToList();
        var phases = lines.Skip(1).Select(l => l.Split('\t')[7]).ToList();

        Assert.That(phases, Is.EqualTo(new[] { "0", "0" }),
            "Each transcript's first CDS segment is phase 0; accumulation must not cross GeneId boundaries");
    }

    /// <summary>
    /// M20: ToGff3 outputs phase "." for non-CDS features.
    /// Source: GFF3 Spec v1.26 — Column 8 phase: 0/1/2 for CDS, "." for others.
    /// NOTE 4: "The phase is REQUIRED for all CDS features."
    /// </summary>
    [Test]
    public void ToGff3_NonCdsFeature_PhaseIsDot()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new(
                GeneId: "gene1",
                Start: 99,
                End: 500,
                Strand: '+',
                Type: "gene",  // non-CDS
                Product: "test",
                Attributes: new Dictionary<string, string>())
        };

        var lines = GenomeAnnotator.ToGff3(annotations, "chr1").ToList();
        var fields = lines[1].Split('\t');

        Assert.That(fields[7], Is.EqualTo("."),
            "Phase must be '.' for non-CDS features per GFF3 Spec v1.26");
    }

    /// <summary>
    /// M21: ToGff3 outputs phase "0" for CDS features.
    /// Source: GFF3 Spec v1.26 — Column 8 phase: 0/1/2 for CDS.
    /// </summary>
    [Test]
    public void ToGff3_CdsFeature_PhaseIsZero()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new(
                GeneId: "cds1",
                Start: 99,
                End: 500,
                Strand: '+',
                Type: "CDS",
                Product: "test protein",
                Attributes: new Dictionary<string, string>())
        };

        var lines = GenomeAnnotator.ToGff3(annotations, "chr1").ToList();
        var fields = lines[1].Split('\t');

        Assert.That(fields[7], Is.EqualTo("0"),
            "Phase must be '0' for CDS features per GFF3 Spec v1.26");
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
            Assert.That(f.Attributes["product"], Is.EqualTo("test protein"), "Product preserved");
            Assert.That(f.Attributes["Note"], Is.EqualTo("important"), "Extra attributes preserved");
        });
    }

    /// <summary>
    /// M22: ToGff3 encodes control characters (tab, newline, CR, percent) per GFF3 Spec.
    /// Source: GFF3 Spec v1.26 — "Literal use of tab, newline, carriage return, the percent (%) sign,
    /// and control characters must be encoded using RFC 3986 Percent-Encoding."
    /// </summary>
    [Test]
    public void ToGff3_EncodesControlCharacters()
    {
        var annotations = new List<GenomeAnnotator.GeneAnnotation>
        {
            new(
                GeneId: "g1",
                Start: 0,
                End: 100,
                Strand: '+',
                Type: "gene",
                Product: "a\tb\nc\r%d",
                Attributes: new Dictionary<string, string>())
        };

        var lines = GenomeAnnotator.ToGff3(annotations).ToList();
        var dataLine = lines[1];

        Assert.That(dataLine, Does.Contain("product=a%09b%0Ac%0D%25d"),
            "Tab, newline, CR, percent must be encoded per GFF3 Spec");
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

    #region B11 review (2026-11-25) — spec conformance vs gt gff3validator / gffutils / BCBio (F8–F10)

    /// <summary>
    /// The Sequence Ontology GFF3 Spec v1.26 "Canonical Gene" (EDEN, lines 2–24), tab-delimited.
    /// It is valid per GenomeTools 1.6.5 <c>gt gff3validator</c>.
    /// </summary>
    private static readonly string[] CanonicalGene =
    {
        "##gff-version 3.1.26",
        "##sequence-region ctg123 1 1497228",
        "ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN",
        "ctg123\t.\tTF_binding_site\t1000\t1012\t.\t+\t.\tID=tfbs00001;Parent=gene00001",
        "ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001;Name=EDEN.1",
        "ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00002;Parent=gene00001;Name=EDEN.2",
        "ctg123\t.\tmRNA\t1300\t9000\t.\t+\t.\tID=mRNA00003;Parent=gene00001;Name=EDEN.3",
        "ctg123\t.\texon\t1300\t1500\t.\t+\t.\tID=exon00001;Parent=mRNA00003",
        "ctg123\t.\texon\t1050\t1500\t.\t+\t.\tID=exon00002;Parent=mRNA00001,mRNA00002",
        "ctg123\t.\texon\t3000\t3902\t.\t+\t.\tID=exon00003;Parent=mRNA00001,mRNA00003",
        "ctg123\t.\texon\t5000\t5500\t.\t+\t.\tID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003",
        "ctg123\t.\texon\t7000\t9000\t.\t+\t.\tID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003",
        "ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1",
        "ctg123\t.\tCDS\t3000\t3902\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1",
        "ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1",
        "ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1",
        "ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00002;Parent=mRNA00002;Name=edenprotein.2",
        "ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00002;Parent=mRNA00002;Name=edenprotein.2",
        "ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00002;Parent=mRNA00002;Name=edenprotein.2",
        "ctg123\t.\tCDS\t3301\t3902\t.\t+\t0\tID=cds00003;Parent=mRNA00003;Name=edenprotein.3",
        "ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00003;Parent=mRNA00003;Name=edenprotein.3",
        "ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00003;Parent=mRNA00003;Name=edenprotein.3",
        "ctg123\t.\tCDS\t3391\t3902\t.\t+\t0\tID=cds00004;Parent=mRNA00003;Name=edenprotein.4",
        "ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00004;Parent=mRNA00003;Name=edenprotein.4",
        "ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00004;Parent=mRNA00003;Name=edenprotein.4",
    };

    private const string ValidLine = "ctg1\t.\tgene\t1\t100\t.\t+\t.\tName=ok";

    private static GenomeAnnotator.GeneAnnotation Ann(
        int start = 9, int end = 20, char strand = '+', string type = "gene", string id = "g1",
        string product = "p", Dictionary<string, string>? attributes = null,
        string source = ".", IReadOnlyDictionary<string, IReadOnlyList<string>>? values = null) =>
        new(id, start, end, strand, type, product, attributes ?? new Dictionary<string, string>(),
            source, null, values);

    private static string[] Row(GenomeAnnotator.GeneAnnotation ann, string seqId = "chr1") =>
        GenomeAnnotator.ToGff3(new[] { ann }, seqId).Skip(1).Single().Split('\t');

    /// <summary>
    /// Spec canonical gene parsed: 23 features; multiple parentage ("using the comma symbol to
    /// indicate multiple parentage") is a value list, identical to gffutils 0.14 / BCBio-GFF
    /// (exon00002 → [mRNA00001, mRNA00002]); seqid/source columns kept.
    /// </summary>
    [Test]
    public void ParseGff3_SpecCanonicalGene_MultipleParentageIsValueList()
    {
        var f = GenomeAnnotator.ParseGff3(CanonicalGene).ToList();

        Assert.That(f, Has.Count.EqualTo(23));
        var exon2 = f.Single(x => x.FeatureId == "exon00002");
        var exon4 = f.Single(x => x.FeatureId == "exon00004");
        Assert.Multiple(() =>
        {
            Assert.That(exon2.AttributeValues!["Parent"], Is.EqualTo(new[] { "mRNA00001", "mRNA00002" }));
            Assert.That(exon4.AttributeValues!["Parent"], Is.EqualTo(new[] { "mRNA00001", "mRNA00002", "mRNA00003" }));
            Assert.That(exon2.Attributes["Parent"], Is.EqualTo("mRNA00001,mRNA00002"));
            Assert.That(f.All(x => x.Seqid == "ctg123" && x.Source == "."), Is.True);
            Assert.That(f.Where(x => x.FeatureId == "cds00003").Select(x => x.Phase), Is.EqualTo(new int?[] { 0, 1, 1 }));
        });
    }

    /// <summary>
    /// F8: an escaped comma belongs to the value, an unescaped comma separates values (GFF3 Spec
    /// column 9). Reference gffutils 0.14: <c>Note=a%2Cb</c> → ['a,b'], <c>Note=a,b</c> → ['a','b'].
    /// Before the fix both lines parsed to the same single string "a,b".
    /// </summary>
    [Test]
    public void ParseGff3_EscapedVersusSeparatorComma_Distinguished()
    {
        var f = GenomeAnnotator.ParseGff3(new[]
        {
            "ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;Note=a%2Cb",
            "ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g2;Note=a,b;Dbxref=EMBL:AA816246,NCBI_gi:10727410",
        }).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(f[0].AttributeValues!["Note"], Is.EqualTo(new[] { "a,b" }));
            Assert.That(f[1].AttributeValues!["Note"], Is.EqualTo(new[] { "a", "b" }));
            Assert.That(f[1].AttributeValues!["Dbxref"], Is.EqualTo(new[] { "EMBL:AA816246", "NCBI_gi:10727410" }));
        });
    }

    /// <summary>
    /// F8: columns 1–3 are kept and percent-decoded; a tag may be escaped ("URL escaping rules are
    /// used for tags or values containing ,=;").
    /// </summary>
    [Test]
    public void ParseGff3_SeqidSourceTypeAndTag_KeptAndDecoded()
    {
        var f = GenomeAnnotator.ParseGff3(new[] { "ctg%3B2\tsrc%09x\tgene\t1\t100\t.\t+\t.\tID=g3;na%3Dme=v" }).Single();

        Assert.Multiple(() =>
        {
            Assert.That(f.Seqid, Is.EqualTo("ctg;2"));
            Assert.That(f.Source, Is.EqualTo("src\tx"));
            Assert.That(f.Attributes["na=me"], Is.EqualTo("v"));
        });
    }

    /// <summary>
    /// F9: "unescaped spaces are allowed within fields" — value text is kept verbatim
    /// (gt gff3 and gffutils keep " x y "; the old parser trimmed it to "x y"); whitespace around a
    /// tag after "; " is not part of the tag (gt, gffutils, BCBio all read "Name").
    /// </summary>
    [Test]
    public void ParseGff3_ValueSpacesKept_TagWhitespaceDropped()
    {
        var f = GenomeAnnotator.ParseGff3(new[] { "ctg1\t.\tgene\t1\t99\t.\t+\t.\tID=g1; Name=x;Note= x y " }).Single();

        Assert.Multiple(() =>
        {
            Assert.That(f.Attributes["Name"], Is.EqualTo("x"));
            Assert.That(f.Attributes["Note"], Is.EqualTo(" x y "));
        });
    }

    /// <summary>
    /// F9: lines that violate the spec's column rules are skipped (one bad line never aborts the
    /// parse and does not consume a feature_n number). Each case is rejected by GenomeTools 1.6.5
    /// <c>gt gff3validator</c> with the quoted error; the old parser accepted all of them.
    /// </summary>
    [TestCase("ctg1\t.\tgene\t1\t100\t.\tZ\t.\tID=x", TestName = "strand 'Z' not a valid character from the set '+-.?'")]
    [TestCase("ctg1\t.\tgene\t1\t100\t.\t\t.\tID=x", TestName = "strand '' not one character long")]
    [TestCase("ctg1\t.\tgene\t1\t100\t.\t+-\t.\tID=x", TestName = "strand '+-' not one character long")]
    [TestCase("ctg1\t.\tCDS\t1\t99\t.\t+\t3\tID=x", TestName = "phase '3' not a valid character from the set '012.'")]
    [TestCase("ctg1\t.\tCDS\t1\t99\t.\t+\t-1\tID=x", TestName = "phase '-1' not one character long")]
    [TestCase("ctg1\t.\tgene\t0\t99\t.\t+\t.\tID=x", TestName = "illegal feature start 0 (GFF3 files are 1-based)")]
    [TestCase("ctg1\t.\tgene\t100\t99\t.\t+\t.\tID=x", TestName = "start '100' is larger then end '99'")]
    [TestCase("ctg1\t.\tgene\t1\t99\t.\t+\t.\tID=x\tID=evil", TestName = "line does not contain 9 tab separated fields (10)")]
    [TestCase("ctg1\t.\tgene\t1\t99\t.\t+\t.", TestName = "line does not contain 9 tab separated fields (8)")]
    [TestCase("ctg1\t.\tgene\t1\t99\t.\t+\t.\tID=x;Note=a;Note=b", TestName = "more than one Note attribute")]
    [TestCase("ctg1\t.\tgene\t1\t99\t.\t+\t.\tID=x;extra", TestName = "token does not contain exactly one '='")]
    [TestCase("ctg1\t.\tgene\t1\t99\t.\t+\t.\tID=x;product=", TestName = "attribute has no value")]
    public void ParseGff3_SpecInvalidLine_Skipped(string invalidLine)
    {
        var f = GenomeAnnotator.ParseGff3(new[] { invalidLine, ValidLine }).ToList();

        Assert.That(f, Has.Count.EqualTo(1));
        Assert.That(f[0].FeatureId, Is.EqualTo("feature_1"), "the skipped line does not consume a feature number");
    }

    /// <summary>
    /// Spec-valid edge lines are accepted (gt gff3validator: "input is valid GFF3"): a zero-length
    /// feature (start = end), the '?' strand, a phase on a non-CDS feature, an E-value score,
    /// an empty column 9 ".", a trailing ';', and a CDS with phase "." (the spec requires a CDS
    /// phase, but gt accepts the line).
    /// </summary>
    [TestCase("ctg1\t.\tinsertion_site\t100\t100\t.\t?\t.\tID=x")]
    [TestCase("ctg1\t.\tgene\t1\t99\t.\t+\t1\tID=x")]
    [TestCase("ctg1\t.\tcDNA_match\t1050\t1500\t5.8e-42\t+\t.\tID=x")]
    [TestCase("ctg1\t.\tgene\t1\t99\t.\t+\t.\t.")]
    [TestCase("ctg1\t.\tgene\t1\t99\t.\t+\t.\tID=x;")]
    [TestCase("ctg1\t.\tCDS\t1\t99\t.\t+\t.\tID=x")]
    public void ParseGff3_SpecValidEdgeLine_Accepted(string line)
    {
        Assert.That(GenomeAnnotator.ParseGff3(new[] { line }).ToList(), Has.Count.EqualTo(1));
    }

    /// <summary>
    /// "##FASTA … indicates that the annotation portion of the file is at an end"; "a GFF line
    /// that begins with the character &gt; creates an implied ##FASTA directive". Nothing after it
    /// is read as a feature, even a tab-containing line.
    /// </summary>
    [TestCase("##FASTA")]
    [TestCase(">ctg1")]
    public void ParseGff3_FastaSection_EndsAnnotation(string fastaStart)
    {
        var f = GenomeAnnotator.ParseGff3(new[]
        {
            "##gff-version 3", ValidLine, fastaStart, "ACGT", "ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=late",
        }).ToList();

        Assert.That(f.Select(x => x.FeatureId), Is.EqualTo(new[] { "feature_1" }));
    }

    /// <summary>A CRLF line terminator is not part of the last value (lines split on '\n').</summary>
    [Test]
    public void ParseGff3_CrLfLine_CrNotInLastValue()
    {
        var f = GenomeAnnotator.ParseGff3(new[] { "ctg1\t.\tgene\t1\t99\t.\t+\t.\tID=g1;Note=x\r" }).Single();
        Assert.That(f.Attributes["Note"], Is.EqualTo("x"));
    }

    [Test]
    public void ParseGff3_NullLines_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => GenomeAnnotator.ParseGff3(null!));
    }

    /// <summary>
    /// The spec canonical gene, rebuilt as GeneAnnotations (0-based start = start − 1, multiple
    /// parentage as AttributeValues), is written back byte-identical to the spec rows — including
    /// all 13 computed CDS phases (cds00003/cds00004: 0,1,1) — and the output is valid per
    /// gt gff3validator.
    /// </summary>
    [Test]
    public void ToGff3_SpecCanonicalGene_ReproducesSpecRows()
    {
        var features = GenomeAnnotator.ParseGff3(CanonicalGene).ToList();
        var annotations = features.Select(f => new GenomeAnnotator.GeneAnnotation(
            f.FeatureId, f.Start - 1, f.End, f.Strand, f.Type, "",
            f.AttributeValues!.Where(kv => kv.Key != "ID" && kv.Value.Count == 1)
                .ToDictionary(kv => kv.Key, kv => kv.Value[0]),
            f.Source, f.Score,
            f.AttributeValues!.Where(kv => kv.Value.Count > 1).ToDictionary(kv => kv.Key, kv => kv.Value))).ToList();

        var lines = GenomeAnnotator.ToGff3(annotations, "ctg123").ToList();

        Assert.That(lines.Skip(1), Is.EqualTo(CanonicalGene.Skip(2)));
    }

    /// <summary>
    /// Phase is measured from the 5' end, which is the feature end on the − strand, so mirroring the
    /// canonical CDS segments onto the − strand (x → L + 1 − x, L = 1497228) keeps every phase:
    /// cds00003 = 0,1,1 in 5'→3' order.
    /// </summary>
    [Test]
    public void ToGff3_CanonicalCdsMirroredToMinusStrand_SamePhases()
    {
        const int L = 1497228;
        (int s, int e)[] cds3 = { (3301, 3902), (5000, 5500), (7000, 7600) };
        var anns = cds3.Select(x => Ann(start: L + 1 - x.e - 1, end: L + 1 - x.s, strand: '-', type: "CDS", id: "cds00003")).ToList();

        var rows = GenomeAnnotator.ToGff3(anns, "ctg123").Skip(1).Select(l => l.Split('\t')).ToList();

        Assert.That(rows.Select(r => r[7]), Is.EqualTo(new[] { "0", "1", "1" }));
        Assert.That(rows.Select(r => r[3] + "-" + r[4]), Is.EqualTo(new[] { "1493327-1493928", "1491729-1492229", "1489629-1490229" }));
    }

    /// <summary>
    /// F10: column 1 "must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|] … may not
    /// contain unescaped whitespace and must not begin with an unescaped '&gt;'" (old output: raw).
    /// </summary>
    [TestCase("chr 1", "chr%201")]
    [TestCase(">chr1", "%3Echr1")]
    [TestCase("chr1;x", "chr1%3Bx")]
    [TestCase("a.b:c^d*e$f@g!h+i_j?k-l|m", "a.b:c^d*e$f@g!h+i_j?k-l|m")]
    [TestCase("ctgé", "ctg%C3%A9")]
    public void ToGff3_SeqidOutsideAllowedSet_PercentEncoded(string seqId, string expected)
    {
        Assert.That(Row(Ann(), seqId)[0], Is.EqualTo(expected));
        Assert.That(GenomeAnnotator.ParseGff3(GenomeAnnotator.ToGff3(new[] { Ann() }, seqId)).Single().Seqid, Is.EqualTo(seqId));
    }

    /// <summary>
    /// F10: tab/newline/CR/'%'/control characters are encoded in every column — a tab in the
    /// source used to split the row into 10 fields (gt: "does not contain 9 tab separated fields").
    /// Column-9 reserved characters are not encoded outside column 9.
    /// </summary>
    [Test]
    public void ToGff3_SourceAndTypeControlChars_Encoded()
    {
        var row = Row(Ann(source: "my\tsrc;1", type: "ge\nne"));

        Assert.That(row, Has.Length.EqualTo(9));
        Assert.That(row[1], Is.EqualTo("my%09src;1"));
        Assert.That(row[2], Is.EqualTo("ge%0Ane"));
    }

    /// <summary>
    /// F10: tags are encoded like values ("URL escaping rules are used for tags or values
    /// containing ',=;'"); the old output "a=b=x" was rejected by gt ("does not contain exactly one '='").
    /// </summary>
    [Test]
    public void ToGff3_TagWithReservedChars_Encoded()
    {
        var row = Row(Ann(attributes: new Dictionary<string, string> { ["a=b;c"] = "x" }));
        Assert.That(row[8], Is.EqualTo("ID=g1;product=p;a%3Db%3Bc=x"));
    }

    /// <summary>
    /// F10: ID and product come only from GeneId/Product — an "ID" key in Attributes used to emit a
    /// second ID tag (gt: "more than one ID attribute"), and ParseGff3 then read the wrong ID.
    /// </summary>
    [Test]
    public void ToGff3_IdAndProductInAttributes_NotRepeated()
    {
        var ann = Ann(attributes: new Dictionary<string, string> { ["ID"] = "x", ["product"] = "y", ["note"] = "n" });

        Assert.That(Row(ann)[8], Is.EqualTo("ID=g1;product=p;note=n"));
        Assert.That(GenomeAnnotator.ParseGff3(GenomeAnnotator.ToGff3(new[] { ann })).Single().FeatureId, Is.EqualTo("g1"));
    }

    /// <summary>
    /// F10: an empty value is omitted (gt: 'attribute "product=" … has no value'); column 9 with no
    /// tag left is the "." placeholder.
    /// </summary>
    [Test]
    public void ToGff3_EmptyValues_Omitted()
    {
        Assert.That(Row(Ann(product: "", attributes: new Dictionary<string, string> { ["note"] = "n" }))[8], Is.EqualTo("ID=g1;note=n"));
        Assert.That(Row(Ann(id: "", product: ""))[8], Is.EqualTo("."));
    }

    /// <summary>
    /// F10: multi-valued attributes ("Parent=AF2312,AB2812,abc-3") are written with unescaped commas
    /// between separately encoded values, and read back as the same list; a single-valued
    /// Attributes entry keeps its comma escaped.
    /// </summary>
    [Test]
    public void ToGff3_AttributeValues_CommaSeparatedAndRoundTrip()
    {
        var values = new Dictionary<string, IReadOnlyList<string>> { ["Parent"] = new[] { "AF2312", "AB2812", "abc-3" }, ["Alias"] = new[] { "x,y", "z" } };
        var ann = Ann(type: "exon", product: "", attributes: new Dictionary<string, string> { ["Note"] = "a,b", ["Parent"] = "ignored" }, values: values);

        var row = Row(ann);
        var f = GenomeAnnotator.ParseGff3(GenomeAnnotator.ToGff3(new[] { ann })).Single();

        Assert.Multiple(() =>
        {
            Assert.That(row[8], Is.EqualTo("ID=g1;Note=a%2Cb;Parent=AF2312,AB2812,abc-3;Alias=x%2Cy,z"));
            Assert.That(f.AttributeValues!["Parent"], Is.EqualTo(new[] { "AF2312", "AB2812", "abc-3" }));
            Assert.That(f.AttributeValues!["Alias"], Is.EqualTo(new[] { "x,y", "z" }));
            Assert.That(f.AttributeValues!["Note"], Is.EqualTo(new[] { "a,b" }));
        });
    }

    /// <summary>
    /// F10: "For zero-length features, such as insertion sites, start equals end and the implied site
    /// is to the right of the indicated base": the empty half-open interval [s, s) is written s..s
    /// (it used to be written start = s + 1 &gt; end = s, rejected by gt).
    /// </summary>
    [Test]
    public void ToGff3_EmptyInterval_WrittenAsZeroLengthFeature()
    {
        var row = Row(Ann(start: 5, end: 5, type: "insertion_site"));
        Assert.That((row[3], row[4]), Is.EqualTo(("5", "5")));
    }

    /// <summary>
    /// F10: annotations without a GFF3 representation are rejected instead of written as invalid
    /// lines (gt: strand 'x' not a valid character from the set '+-.?'; start larger than end).
    /// </summary>
    [TestCase(9, 20, 'x', "gene")]
    [TestCase(9, 20, '\0', "gene")]
    [TestCase(-1, 20, '+', "gene")]
    [TestCase(20, 9, '+', "gene")]
    [TestCase(0, 0, '+', "gene")]
    [TestCase(9, 20, '+', "")]
    public void ToGff3_UnrepresentableAnnotation_Throws(int start, int end, char strand, string type)
    {
        Assert.Throws<ArgumentException>(() => GenomeAnnotator.ToGff3(new[] { Ann(start, end, strand, type) }).ToList());
    }

    [Test]
    public void ToGff3_NullAnnotationsOrEmptySeqId_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => GenomeAnnotator.ToGff3(null!));
        Assert.Throws<ArgumentException>(() => GenomeAnnotator.ToGff3(new[] { Ann() }, ""));
    }

    /// <summary>The '?' and '.' strands are valid GFF3 column-7 symbols and are written verbatim.</summary>
    [TestCase('?')]
    [TestCase('.')]
    public void ToGff3_UnknownAndUnstrandedSymbols_Written(char strand)
    {
        Assert.That(Row(Ann(strand: strand))[6], Is.EqualTo(strand.ToString()));
    }

    #endregion
}
