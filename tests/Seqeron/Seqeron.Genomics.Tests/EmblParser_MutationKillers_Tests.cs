using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.IO;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Targeted mutation-killing tests for EmblParser.cs (checklist 04 row 70, PARSE-EMBL-001).
/// The canonical suite left the ID-line decomposition, line-type grouping, multi-block reference
/// assembly, multi-line FT location/qualifier continuation and the INSDC quote-unescaping
/// under-pinned. These pin the exact EBI EMBL flat-file behaviour so the boundary/logical mutants
/// diverge.
/// </summary>
[TestFixture]
public class EmblParser_MutationKillers_Tests
{
    private static readonly string Record = string.Join("\n", new[]
    {
        "ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.",
        "XX",
        "AC   X56734; S46826;",
        "XX",
        "DE   Trifolium repens mRNA for beta-glucosidase",
        "XX",
        "KW   beta-glucosidase.",
        "XX",
        "OS   Trifolium repens (white clover)",
        "OC   Eukaryota; Viridiplantae; Streptophyta.",
        "XX",
        "RN   [1]",
        "RP   1-1859",
        "RA   Oxtoby E., Dunn M.A.;",
        "RT   \"Nucleotide sequence\";",
        "RL   Plant Mol. Biol. 17:209-219(1991).",
        "XX",
        "RN   [2]",
        "RA   Brown,T.;",
        "XX",
        "FH   Key             Location/Qualifiers",
        "FT   gene            14..898",
        "FT                   /gene=\"bgl\"",
        "FT   CDS             join(14..200,",
        "FT                   300..898)",
        "FT                   /gene=\"bgl\"",
        "FT                   /product=\"beta-",
        "FT                   glucosidase\"",
        "FT                   /number=12",
        "FT                   /note=\"a \"\"b\"\" c\"",
        "FT                   /experiment=\"\"",
        "XX",
        "SQ   Sequence 1859 BP;",
        "     aaacaaacca aatatggatt                                              20",
        "//"
    });

    // ── ID line decomposition ───────────────────────────────────────────────────────────

    [Test]
    public void Parse_IdLine_DecomposesAllTokens()
    {
        var rec = EmblParser.Parse(Record).Single();

        rec.Accession.Should().Be("X56734");
        rec.SequenceVersion.Should().Be("1");
        rec.Topology.Should().Be("linear");
        rec.MoleculeType.Should().Be("mRNA");
        rec.DataClass.Should().Be("STD");
        rec.TaxonomicDivision.Should().Be("PLN");
        rec.SequenceLength.Should().Be(1859);
    }

    // ── Descriptive line groups ─────────────────────────────────────────────────────────

    [Test]
    public void Parse_DescriptionKeywordsOrganismClassification_ParsedExactly()
    {
        var rec = EmblParser.Parse(Record).Single();

        rec.Description.Should().Be("Trifolium repens mRNA for beta-glucosidase");
        rec.Keywords.Should().Equal("beta-glucosidase");
        rec.Organism.Should().Be("Trifolium repens (white clover)");
        rec.OrganismClassification.Should().Equal("Eukaryota", "Viridiplantae", "Streptophyta");
    }

    // ── Sequence (SQ … //) — letters only, uppercased, digits/whitespace dropped ───────

    [Test]
    public void Parse_Sequence_KeepsOnlyUppercasedLetters()
    {
        EmblParser.Parse(Record).Single().Sequence.Should().Be("AAACAAACCAAATATGGATT");
    }

    // ── Reference assembly across two RN blocks ────────────────────────────────────────

    [Test]
    public void Parse_References_AssemblesBothBlocks()
    {
        var refs = EmblParser.Parse(Record).Single().References;

        // Both references must be retained (currentRefNum > 0 on save — kills the `< 0` mutants on
        // the per-RN flush and the final save).
        refs.Should().HaveCount(2);
        refs[0].Number.Should().Be(1);
        refs[0].Positions.Should().Be("1-1859");
        refs[0].Authors.Should().Be("Oxtoby E., Dunn M.A.");
        refs[0].Title.Should().Be("Nucleotide sequence");
        refs[0].Journal.Should().Be("Plant Mol. Biol. 17:209-219(1991).");

        refs[1].Number.Should().Be(2);
        refs[1].Authors.Should().Be("Brown,T.");
    }

    // ── Accession/version fall back to dedicated AC / SV lines when absent from ID ─────

    [Test]
    public void Parse_AccessionAndVersionFallBackToAcAndSvLines()
    {
        // ID line carries no accession token and no "SV n" ⇒ the parser must source the
        // accession from the AC line and the version from the SV line (covers the two fallback
        // branches the canonical sample never reached).
        var content = string.Join("\n", new[]
        {
            "ID   ; ; linear; mRNA; STD; PLN; 100 BP.",
            "XX",
            "AC   Y99999; Z11111;",
            "XX",
            "SV   Y99999.2",
            "XX",
            "SQ   Sequence 100 BP;",
            "     acgtacgtac",
            "//"
        });

        var rec = EmblParser.Parse(content).Single();
        rec.Accession.Should().Be("Y99999");        // primary accession from AC line
        rec.SequenceVersion.Should().Be("Y99999.2"); // version from SV line
    }

    // ── FT feature: multi-line location & qualifiers, INSDC quote handling ─────────────

    [Test]
    public void Parse_Feature_HandlesMultiLineLocationAndQualifiers()
    {
        var rec = EmblParser.Parse(Record).Single();

        // Two features: the gene (flushed mid-stream when the CDS line starts) and the CDS.
        rec.Features.Should().HaveCount(2);
        rec.Features[0].Key.Should().Be("gene");
        rec.Features[0].Location.RawLocation.Should().Be("14..898");
        rec.Features[0].Qualifiers["gene"].Should().Be("bgl");

        var cds = rec.Features[1];
        cds.Key.Should().Be("CDS");
        // Location continued on the next FT line (kills the location-continuation guards):
        cds.Location.RawLocation.Should().Be("join(14..200,300..898)");

        cds.Qualifiers["gene"].Should().Be("bgl");
        // Quoted value split across two FT lines (kills the qualifier-continuation guards):
        cds.Qualifiers["product"].Should().Be("beta-glucosidase");
        // Unquoted value must NOT be quote-stripped (kills the `value.Length >= 2 || …` logical mutant):
        cds.Qualifiers["number"].Should().Be("12");
        // INSDC: embedded quote encoded as doubled quote, collapsed on parse:
        cds.Qualifiers["note"].Should().Be("a \"b\" c");
        // Empty quoted value unquotes to "" (kills the `value.Length > 2` boundary):
        cds.Qualifiers["experiment"].Should().Be("");
    }
}
