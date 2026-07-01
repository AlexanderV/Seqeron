using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.IO;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// Targeted mutation-killing tests for GenBankParser.cs (checklist 04 row 69, PARSE-GENBANK-001).
/// The canonical suite left the LOCUS line decomposition, section extraction, multi-block
/// reference assembly, multi-line feature qualifiers, the standard-field partition and the CDS
/// codon-translation table under-pinned. These pin the exact NCBI/INSDC flat-file behaviour and
/// the standard genetic code so the boundary/logical/switch mutants diverge.
/// </summary>
[TestFixture]
public class GenBankParserMutationTests
{
    private static readonly string Record = string.Join("\n", new[]
    {
        "LOCUS       SCU49845     5028 bp    DNA     linear   PLN 21-JUN-1999",
        "DEFINITION  Saccharomyces cerevisiae TCP1-beta gene.",
        "ACCESSION   U49845",
        "VERSION     U49845.1",
        "KEYWORDS    beta; TCP1.",
        "SOURCE      Saccharomyces cerevisiae (baker's yeast)",
        "  ORGANISM  Saccharomyces cerevisiae",
        "            Eukaryota; Fungi; Ascomycota;",
        "            Saccharomycetales.",
        "REFERENCE   1  (bases 1 to 5028)",
        "  AUTHORS   Roemer,T.",
        "  TITLE     Selection of axial growth sites in yeast",
        "  JOURNAL   Cell 79, 1069-1080 (1994)",
        "  PUBMED    7954430",
        "REFERENCE   2",
        "  AUTHORS   Smith,J.",
        "COMMENT     This is a test comment.",
        "DBLINK      BioProject: PRJNA000",
        "FEATURES             Location/Qualifiers",
        "     gene            1..206",
        "                     /gene=\"TCP1\"",
        "     CDS             1..12",
        "                     /gene=\"TCP1\"",
        "                     /product=\"TCP1-beta\"",
        "                     /function=cell wall",
        "                     biogenesis",
        "                     /note=\"\"",
        "NID         g1293613",
        "ORIGIN",
        "        1 atggtgcacc tgactcctga ggagaagtct",
        "//"
    });

    // ── A complete record decomposes into every field exactly ──────────────────────────

    [Test]
    public void Parse_FullRecord_DecomposesEveryFieldExactly()
    {
        var rec = GenBankParser.Parse(Record).Single();

        rec.Locus.Should().Be("SCU49845");
        rec.SequenceLength.Should().Be(5028);
        rec.MoleculeType.Should().Be("DNA");
        rec.Topology.Should().Be("linear");
        rec.Division.Should().Be("PLN");
        rec.Date.Should().Be(new System.DateTime(1999, 6, 21));

        rec.Definition.Should().Be("Saccharomyces cerevisiae TCP1-beta gene.");
        rec.Accession.Should().Be("U49845");
        rec.Version.Should().Be("U49845.1");
        rec.Keywords.Should().Equal("beta", "TCP1");
        rec.Organism.Should().Be("Saccharomyces cerevisiae");
        // Two taxonomy continuation lines joined by a single space (kills the `taxonomy.Length > 0`
        // separator guard, which would concatenate them with no space).
        rec.Taxonomy.Should().Be("Eukaryota; Fungi; Ascomycota; Saccharomycetales");

        rec.Sequence.Should().Be("ATGGTGCACCTGACTCCTGAGGAGAAGTCT");
    }

    [Test]
    public void Parse_References_AssemblesBothBlocksIncludingTitlelessOne()
    {
        var refs = GenBankParser.Parse(Record).Single().References;

        // The second reference has only AUTHORS (no TITLE); it must still be retained because
        // refNum > 0 (kills the `refNum > 0 || title` → `&&` and `refNum < 0` mutants).
        refs.Should().HaveCount(2);
        refs[0].Number.Should().Be(1);
        refs[0].BaseFrom.Should().Be(1);
        refs[0].BaseTo.Should().Be(5028);
        refs[0].Authors.Should().Be("Roemer,T.");
        refs[0].Title.Should().Be("Selection of axial growth sites in yeast");
        refs[0].Journal.Should().Be("Cell 79, 1069-1080 (1994)");
        refs[0].PubMed.Should().Be("7954430");

        refs[1].Number.Should().Be(2);
        refs[1].Authors.Should().Be("Smith,J.");
        refs[1].Title.Should().Be("");
    }

    [Test]
    public void Parse_Features_ParsesKeysLocationsAndMultiLineAndEmptyQualifiers()
    {
        var rec = GenBankParser.Parse(Record).Single();

        rec.Features.Should().HaveCount(2);
        rec.Features[0].Key.Should().Be("gene");
        rec.Features[0].Location.RawLocation.Should().Be("1..206");
        rec.Features[0].Qualifiers["gene"].Should().Be("TCP1");

        var cds = rec.Features[1];
        cds.Key.Should().Be("CDS");
        cds.Location.RawLocation.Should().Be("1..12");
        cds.Qualifiers["gene"].Should().Be("TCP1");
        cds.Qualifiers["product"].Should().Be("TCP1-beta");
        // Multi-line qualifier value concatenated across the continuation line:
        cds.Qualifiers["function"].Should().Be("cell wall biogenesis");
        // Empty quoted value "" unquotes to the empty string (kills `value.Length > 2` boundary).
        cds.Qualifiers["note"].Should().Be("");
    }

    [Test]
    public void Parse_AdditionalFields_HoldOnlyNonStandardSections()
    {
        var rec = GenBankParser.Parse(Record).Single();

        // Standard sections route to typed fields; only the non-standard "NID" lands in
        // AdditionalFields. This pins the IsStandardField partition.
        rec.AdditionalFields.Keys.Should().BeEquivalentTo(new[] { "NID" });
        rec.AdditionalFields["NID"].Should().Be("g1293613");
    }

    // ── TranslateCDS: full standard genetic code (kills the codon switch) ──────────────

    [Test]
    public void TranslateCDS_TranslatesAll64CodonsPerStandardGeneticCode()
    {
        const string bases = "ACGT";
        var dna = new StringBuilder();
        foreach (var b1 in bases)
            foreach (var b2 in bases)
                foreach (var b3 in bases)
                    dna.Append(b1).Append(b2).Append(b3);

        const string expected =
            "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

        GenBankParser.TranslateCDS(BuildCdsRecord(dna.ToString()), BuildCds(dna.ToString())).Should().Be(expected);
    }

    [Test]
    public void TranslateCDS_TrailingPartialCodonIsNotTranslated()
    {
        // 5-base CDS ⇒ one full codon (ATG→M) plus a 2-base remainder that the loop bound
        // `i + 2 < dna.Length` must exclude (kills `<` → `<=`, which would read past the string).
        GenBankParser.TranslateCDS(BuildCdsRecord("ATGCC"), BuildCds("ATGCC")).Should().Be("M");
    }

    [Test]
    public void TranslateCDS_PrefersExistingTranslationQualifier()
    {
        var cds = new GenBankParser.Feature("CDS", GenBankParser.ParseLocation("1..3"),
            new Dictionary<string, string> { ["translation"] = "MKV" });
        var record = new GenBankParser.GenBankRecord(
            "L", 3, "DNA", "linear", "PLN", null, "", "", "",
            System.Array.Empty<string>(), "", "",
            System.Array.Empty<GenBankParser.Reference>(),
            new[] { cds }, "ATG", new Dictionary<string, string>());

        GenBankParser.TranslateCDS(record, cds).Should().Be("MKV");
    }

    private static GenBankParser.Feature BuildCds(string seq)
        => new("CDS", GenBankParser.ParseLocation($"1..{seq.Length}"), new Dictionary<string, string>());

    private static GenBankParser.GenBankRecord BuildCdsRecord(string seq)
        => new("L", seq.Length, "DNA", "linear", "PLN", null, "", "", "",
            System.Array.Empty<string>(), "", "",
            System.Array.Empty<GenBankParser.Reference>(),
            new[] { BuildCds(seq) }, seq, new Dictionary<string, string>());
}
