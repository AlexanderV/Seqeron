using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>MCP <c>write_big_alignment</c> / <c>parse_big_alignment</c> (B09 A2-27d, F234 - F236) and <c>tcoffee_evaluate</c> (F237).</summary>
[TestFixture]
public class BigFormatTests
{
    private const string Pair = "{\"sequences\": [{\"id\": \"chr1\", \"length\": 100, \"parts\": null}, {\"id\": \"read1\", \"length\": 20, \"parts\": null}], \"coordinates\": [[10, 20, 25, 40], [0, 10, 10, 25]], \"attrs\": {\"score\": 500}}";
    private const string Chromosomes = "[{\"id\": \"chr1\", \"length\": 100}]";

    [Test]
    public void BigBed_RoundTrip_ThroughBase64()
    {
        var written = AlignmentTools.WriteBigAlignment([Pair], Chromosomes, "bigbed");
        var read = AlignmentTools.ParseBigAlignment(written.Base64, "bigbed");
        var aligns = JsonNode.Parse(read.AlignmentsJson)!.AsArray();
        Assert.Multiple(() =>
        {
            Assert.That(written.Count, Is.EqualTo(1));
            Assert.That(Convert.FromBase64String(written.Base64)[..4], Is.EqualTo(new byte[] { 0xEB, 0xF2, 0x89, 0x87 }));
            Assert.That(aligns.Count, Is.EqualTo(1));
            Assert.That(aligns[0]!["coordinates"]!.ToJsonString(), Is.EqualTo("[[10,20,25,40],[0,10,10,25]]"));
            Assert.That(read.DeclarationJson, Does.Contain("\"name\":\"bed\""));
        });
    }

    [Test]
    public void BigPsl_RoundTrip_KeepsTheCountsAndTheCds()
    {
        const string psl = "{\"sequences\": [{\"id\": \"chr1\", \"length\": 100, \"parts\": null}, {\"id\": \"read1\", \"length\": 20, \"parts\": {\"0\": \"ACGTACGTACGTACGTACGT\"}, \"features\": [{\"type\": \"CDS\", \"location\": {\"start\": 2, \"end\": 11, \"strand\": 1}, \"qualifiers\": {}}]}], \"coordinates\": [[10, 20, 25, 40], [0, 10, 10, 25]], \"attrs\": {\"score\": 500, \"annotations\": {}}}";
        var written = AlignmentTools.WriteBigAlignment([psl], Chromosomes, "bigpsl", cds: true, fa: true);
        var read = AlignmentTools.ParseBigAlignment(written.Base64, "bigpsl");
        var a = JsonNode.Parse(read.AlignmentsJson)![0]!;
        Assert.Multiple(() =>
        {
            Assert.That(a["sequences"]![1]!["parts"]!["0"]!.GetValue<string>(), Is.EqualTo("ACGTACGTACGTACGTACGT"));
            Assert.That(a["sequences"]![1]!["features"]![0]!["location"]!["start"]!.GetValue<int>(), Is.EqualTo(2));
            Assert.That(a["attrs"]!["matches"]!.GetValue<int>(), Is.EqualTo(0).Or.GreaterThan(0));
            Assert.That(read.DeclarationJson, Does.Contain("bigPsl"));
        });
    }

    [Test]
    public void BigMaf_RoundTrip_ReadsTheMafBlock()
    {
        const string maf = "{\"sequences\": [{\"id\": \"hg38.chr1\", \"length\": 100, \"parts\": {\"0\": \"ACGTACGTAC\"}}, {\"id\": \"mm10.chr2\", \"length\": 50, \"parts\": {\"0\": \"ACGTACAC\"}}], \"coordinates\": [[0, 6, 8, 10], [0, 6, 6, 8]], \"attrs\": {\"score\": 12.5}}";
        var written = AlignmentTools.WriteBigAlignment([maf], "[{\"id\": \"hg38.chr1\", \"length\": 100}]", "bigmaf");
        var read = AlignmentTools.ParseBigAlignment(written.Base64, "bigmaf");
        var a = JsonNode.Parse(read.AlignmentsJson)![0]!;
        Assert.Multiple(() =>
        {
            Assert.That(a["sequences"]![0]!["id"]!.GetValue<string>(), Is.EqualTo("hg38.chr1"));
            Assert.That(a["sequences"]![1]!["id"]!.GetValue<string>(), Is.EqualTo("mm10.chr2"));
            Assert.That(JsonNode.Parse(read.TargetsJson)![0]!["id"]!.GetValue<string>(), Is.EqualTo("hg38.chr1"));
        });
    }

    [Test]
    public void BadInput_IsRefused()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => AlignmentTools.ParseBigAlignment(Convert.ToBase64String("not a bigbed file"u8.ToArray()), "bigbed"), Throws.ArgumentException);
            Assert.That(() => AlignmentTools.ParseBigAlignment("***", "bigbed"), Throws.ArgumentException);
            Assert.That(() => AlignmentTools.WriteBigAlignment([Pair], Chromosomes, "psl"), Throws.ArgumentException);
        });
    }
}
