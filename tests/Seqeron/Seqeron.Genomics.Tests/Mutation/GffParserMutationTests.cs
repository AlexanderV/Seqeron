using static Seqeron.Genomics.IO.GffParser;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// Targeted mutation-killing tests for GffParser.cs (checklist 04 row 68, PARSE-GFF-001).
/// The canonical suite left the GFF3 1-based coordinate handling, the gene-name attribute
/// precedence (gene_name → Name → gene_id → ID), region overlap bounds and overlap-merge
/// adjacency under-pinned. These pin the exact spec behaviour so the boundary/null-coalescing
/// mutants diverge.
/// </summary>
[TestFixture]
public class GffParserMutationTests
{
    private const string Gff3 =
        "##gff-version 3\n" +
        "chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tID=g1;Name=NAME1;gene_name=GN1\n" +
        "chr1\tHAVANA\tgene\t3000\t4000\t12.5\t-\t.\tID=g2;Name=NAME2\n" +
        "chr1\tHAVANA\tgene\t5000\t6000\t.\t+\t.\tID=g3;gene_id=GI3\n" +
        "chr1\tHAVANA\tgene\t7000\t8000\t.\t+\t.\tID=g4\n";

    // ── ParseLine: every column parsed exactly ─────────────────────────────────────────

    [Test]
    public void Parse_Gff3_ParsesAllColumnsExactly()
    {
        var records = Parse(Gff3).ToList();

        records.Should().HaveCount(4);
        var g1 = records[0];
        g1.Seqid.Should().Be("chr1");
        g1.Source.Should().Be("HAVANA");
        g1.Type.Should().Be("gene");
        g1.Start.Should().Be(1000);
        g1.End.Should().Be(2000);
        g1.Score.Should().BeNull();          // "." ⇒ null score
        g1.Strand.Should().Be('+');
        g1.Phase.Should().BeNull();

        records[1].Score.Should().Be(12.5);   // numeric score parsed
        records[1].Strand.Should().Be('-');
    }

    // ── GetGeneName: gene_name → Name → gene_id → ID coalescing precedence ─────────────

    [Test]
    public void GetGeneName_FollowsAttributePrecedence()
    {
        var r = Parse(Gff3).ToList();

        GetGeneName(r[0]).Should().Be("GN1");   // gene_name wins over Name
        GetGeneName(r[1]).Should().Be("NAME2"); // no gene_name ⇒ Name
        GetGeneName(r[2]).Should().Be("GI3");   // no gene_name/Name ⇒ gene_id (over ID)
        GetGeneName(r[3]).Should().Be("g4");    // only ID present
    }

    // ── FilterByRegion: inclusive interval-overlap test (Start <= end && End >= start) ──

    [Test]
    public void FilterByRegion_OverlapBoundsAreInclusive()
    {
        var records = new[]
        {
            new GffRecord("chr1", "s", "x", 10, 20, null, '+', null, new Dictionary<string, string> { ["ID"] = "touchEnd" }),
            new GffRecord("chr1", "s", "x", 30, 40, null, '+', null, new Dictionary<string, string> { ["ID"] = "touchStart" }),
            new GffRecord("chr1", "s", "x", 50, 60, null, '+', null, new Dictionary<string, string> { ["ID"] = "outside" }),
        };

        // Query [20,30]: feature ending exactly at 20 overlaps (End >= start boundary),
        // feature starting exactly at 30 overlaps (Start <= end boundary), [50,60] does not.
        FilterByRegion(records, "chr1", 20, 30)
            .Select(r => r.Attributes["ID"])
            .Should().BeEquivalentTo("touchEnd", "touchStart");
    }

    // ── ExtractSequence: GFF is 1-based; minus strand returns reverse complement ────────

    [Test]
    public void ExtractSequence_PlusStrand_UsesOneBasedCoordinates()
    {
        var rec = new GffRecord("chr1", "s", "x", 2, 4, null, '+', null, new Dictionary<string, string>());
        // 1-based [2,4] over "ACGTAC" ⇒ 0-based slice [1,4) = "CGT".
        ExtractSequence(rec, "ACGTAC").Should().Be("CGT");
    }

    [Test]
    public void ExtractSequence_MinusStrand_ReverseComplements()
    {
        var rec = new GffRecord("chr1", "s", "x", 2, 4, null, '-', null, new Dictionary<string, string>());
        // reverse complement of "CGT" = "ACG".
        ExtractSequence(rec, "ACGTAC").Should().Be("ACG");
    }

    // ── MergeOverlapping: features touching within 1 bp are merged (Start <= End + 1) ──

    [Test]
    public void MergeOverlapping_AdjacentFeaturesMergeAcrossOneBpGap()
    {
        var records = new[]
        {
            new GffRecord("chr1", "s", "x", 1, 10, null, '+', null, new Dictionary<string, string>()),
            new GffRecord("chr1", "s", "x", 11, 20, null, '+', null, new Dictionary<string, string>()),
        };

        var merged = MergeOverlapping(records).ToList();

        // next.Start (11) <= current.End (10) + 1 ⇒ merge into a single [1,20] feature.
        merged.Should().HaveCount(1);
        merged[0].Start.Should().Be(1);
        merged[0].End.Should().Be(20);
    }

    [Test]
    public void MergeOverlapping_GapLargerThanOne_KeepsSeparate()
    {
        var records = new[]
        {
            new GffRecord("chr1", "s", "x", 1, 10, null, '+', null, new Dictionary<string, string>()),
            new GffRecord("chr1", "s", "x", 12, 20, null, '+', null, new Dictionary<string, string>()),
        };

        // next.Start (12) > current.End (10) + 1 ⇒ NOT merged (kills `<= End + 1` → `< End + 1`
        // and the `End - 1` arithmetic mutant, which would still merge).
        MergeOverlapping(records).Should().HaveCount(2);
    }

    // ── Round-trip: WriteToStream(GFF3) then Parse recovers fields and attributes ───────

    [Test]
    public void WriteThenParse_Gff3_RoundTripsRecord()
    {
        var rec = new GffRecord("chrX", "src", "exon", 100, 200, 9.0, '-', 0,
            new Dictionary<string, string> { ["ID"] = "e1", ["Parent"] = "t1" });

        using var sw = new StringWriter();
        WriteToStream(sw, new[] { rec }, GffFormat.GFF3);
        var text = sw.ToString();

        text.Should().StartWith("##gff-version 3");   // kills the GFF3 header guard

        var parsed = Parse(text).Single();
        parsed.Seqid.Should().Be("chrX");
        parsed.Type.Should().Be("exon");
        parsed.Start.Should().Be(100);
        parsed.End.Should().Be(200);
        parsed.Strand.Should().Be('-');
        parsed.Phase.Should().Be(0);
        parsed.Attributes["ID"].Should().Be("e1");
        parsed.Attributes["Parent"].Should().Be("t1");
    }
}
