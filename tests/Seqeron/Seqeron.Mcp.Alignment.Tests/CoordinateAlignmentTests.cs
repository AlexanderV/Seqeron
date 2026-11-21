using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>MCP <c>write_coordinate_alignment</c> / <c>parse_coordinate_alignment</c> (B09 A2-27c, F224): text equal to Biopython 1.88 (<c>Bio.Align.write</c> of the same alignment) and the round trip.</summary>
[TestFixture]
public class CoordinateAlignmentTests
{
    private const string Alignment = "{\"sequences\":[{\"id\":\"chr1\",\"length\":20,\"parts\":{\"0\":\"ACGTACGTAAGGTTACGTAC\"}},{\"id\":\"q1\",\"length\":12,\"parts\":{\"0\":\"ACGTACTAAGGT\"}}],\"coordinates\":[[2,8,10,16],[0,6,6,12]],\"attrs\":{\"score\":42}}";

    [Test]
    public void Write_Bed_Chain_Psl_EqualBiopython()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AlignmentTools.WriteCoordinateAlignment([Alignment], "bed", bedN: 6).Text, Is.EqualTo("chr1\t2\t16\tq1\t42\t+\n"));
            Assert.That(AlignmentTools.WriteCoordinateAlignment([Alignment], "chain").Text, Is.EqualTo("chain 42 chr1 20 + 2 16 q1 12 + 0 12\n6\t2\t0\n6\n\n"));
            var psl = AlignmentTools.WriteCoordinateAlignment([Alignment], "psl", header: false);
            Assert.That(psl.Text, Is.EqualTo("0\t12\t0\t0\t0\t0\t1\t2\t+\tq1\t12\t0\t12\tchr1\t20\t2\t16\t2\t6,6,\t0,6,\t2,10,\n"));
            Assert.That(psl.Count, Is.EqualTo(1));
        });
    }

    [Test]
    public void Parse_RoundTripsTheChainBlock()
    {
        var parsed = AlignmentTools.ParseCoordinateAlignment("chain 42 chr1 20 + 2 16 q1 12 + 0 12\n6\t2\t0\n6\n\n", "chain");
        var arr = JsonNode.Parse(parsed.Text)!.AsArray();
        Assert.Multiple(() =>
        {
            Assert.That(arr, Has.Count.EqualTo(1));
            Assert.That(arr[0]!["coordinates"]!.ToJsonString(), Is.EqualTo("[[2,8,10,16],[0,6,6,12]]"));
            Assert.That(arr[0]!["attrs"]!["score"]!.GetValue<double>(), Is.EqualTo(42.0));
        });
    }

    [Test]
    public void Write_And_Parse_Maf_EqualBiopython()
    {
        const string maf = "{\"sequences\":[{\"id\":\"hg38.chr1\",\"length\":12,\"parts\":{\"0\":\"ACGTACGTAAGG\"}},{\"id\":\"mm10.chr2\",\"length\":9,\"parts\":{\"0\":\"ACTACGAAG\"}}],\"coordinates\":[[2,5,6,11],[1,4,4,9]],\"attrs\":{\"score\":42}}";
        string expected = "##maf version=1\na score=42.000000\ns hg38.chr1 2 9 + 12 GTACGTAAG\ns mm10.chr2 1 8 +  9 CTA-CGAAG\n\n";
        var written = AlignmentTools.WriteCoordinateAlignment([maf], "maf");
        var parsed = JsonNode.Parse(AlignmentTools.ParseCoordinateAlignment(written.Text, "maf").Text)!.AsArray();
        Assert.Multiple(() =>
        {
            Assert.That(written.Text, Is.EqualTo(expected));
            Assert.That(parsed[0]!["coordinates"]!.ToJsonString(), Is.EqualTo("[[2,5,6,11],[1,4,4,9]]"));
            Assert.That(parsed[0]!["attrs"]!["score"]!.GetValue<double>(), Is.EqualTo(42.0));
        });
    }

    [Test]
    public void UnknownFormat_IsRefused()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.WriteCoordinateAlignment([Alignment], "sam"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ParseCoordinateAlignment("x", "sam"));
    }
}
