using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-27d, B09 F234): the bigBed writer (<see cref="AlignmentCoordinateFormats.WriteBigBed"/>) writes the bytes of <c>Bio.Align.write(alignments, "bigbed")</c> (the data blocks are compressed with the zlib calls Biopython
/// made, which depend on the zlib build) and the reader (<see cref="AlignmentCoordinateFormats.ParseBigBed"/>) equals <c>Bio.Align.parse</c> (scripts/alignment/oracle_sections_bigbed.py, biopython_coord_oracle.tsv).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class AlignmentCoordinateBigBed_Tests
{
    private static readonly Lazy<List<BiopythonCoordOracle.Row>> Oracle = new(() => BiopythonCoordOracle.Load().Where(r => r.Format == "bigbed").ToList());

    private static BigBedWriteOptions Options(JsonObject o, Dictionary<string, byte[]> compressed)
    {
        IReadOnlyList<BigBedField>? fields = null;
        string name = "bed", comment = "Browser Extensible Data";
        if (o["declaration"] is JsonObject d)
        {
            name = d["name"]!.GetValue<string>();
            comment = d["comment"]!.GetValue<string>();
            fields = [.. d["fields"]!.AsArray().Select(f => new BigBedField(f![0]!.GetValue<string>(), f[1]!.GetValue<string>(), f[2]!.GetValue<string>()))];
        }

        return new BigBedWriteOptions
        {
            BedN = o["bedN"]?.GetValue<int>() ?? 12,
            Compress = o["compress"]?.GetValue<bool>() ?? true,
            ItemsPerSlot = o["itemsPerSlot"]?.GetValue<int>() ?? 512,
            BlockSize = o["blockSize"]?.GetValue<int>() ?? 256,
            ExtraIndex = o["extraIndex"] is JsonArray ei ? [.. ei.Select(x => x!.GetValue<string>())] : [],
            DeclarationFields = fields,
            DeclarationName = name,
            DeclarationComment = comment,
            Compressor = data => compressed.TryGetValue(Convert.ToBase64String(data), out var c) ? c : throw new InvalidOperationException("a block that Biopython did not compress"),
        };
    }

    [Test]
    public void Oracle_Has200WriterAndReaderRows()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Oracle.Value.Count(r => r.Kind == "write" && !r.Data2.StartsWith("{\"error", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(200));
            Assert.That(Oracle.Value.Count(r => r.Kind == "read"), Is.GreaterThanOrEqualTo(200));
        });
    }

    [Test]
    public void Write_EveryRow_IsByteEqualToBiopython()
    {
        var failures = new List<string>();
        foreach (var r in Oracle.Value.Where(r => r.Kind == "write"))
        {
            var input = JsonNode.Parse(r.Data1)!;
            var aligns = input["alignments"]!.AsArray().Select(s => BiopythonCoordOracle.FromSpec(s!)).ToList();
            var targets = input["targets"]!.AsArray().Select(t => AlignmentCoordinateFormats.FromJson(new JsonObject { ["sequences"] = new JsonArray(t!.DeepClone()), ["coordinates"] = new JsonArray(), ["attrs"] = new JsonObject() }).Sequences[0]).ToList();
            var expected = JsonNode.Parse(r.Data2)!.AsObject();
            var compressed = new Dictionary<string, byte[]>();
            foreach (var p in expected["zlib"]?.AsArray() ?? [])
                compressed[p![0]!.GetValue<string>()] = Convert.FromBase64String(p[1]!.GetValue<string>());
            try
            {
                string actual = Convert.ToBase64String(AlignmentCoordinateFormats.WriteBigBed(aligns, targets, Options(r.Options, compressed)));
                if (expected.ContainsKey("error"))
                    failures.Add($"{r.Id} wrote but Biopython refused ({expected["error"]})");
                else if (actual != expected["file"]!.GetValue<string>())
                    failures.Add($"{r.Id} {r.Options.ToJsonString()}: file differs (length {Convert.FromBase64String(actual).Length} / {Convert.FromBase64String(expected["file"]!.GetValue<string>()).Length})");
            }
            catch (ArgumentException)
            {
                if (!expected.ContainsKey("error"))
                    failures.Add($"{r.Id} {r.Options.ToJsonString()} threw but Biopython wrote");
            }
        }

        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(6)) + $"\n({failures.Count})");
    }

    [Test]
    public void Parse_EveryRow_EqualsBioAlignParse()
    {
        var failures = new List<string>();
        foreach (var r in Oracle.Value.Where(r => r.Kind == "read"))
        {
            var expected = JsonNode.Parse(r.Data2)!.AsObject();
            CoordinateAlignmentFile file;
            try
            {
                file = AlignmentCoordinateFormats.ParseBigBed(Convert.FromBase64String(r.Data1));
            }
            catch (ArgumentException)
            {
                if (!expected.ContainsKey("error"))
                    failures.Add($"{r.Id} threw but Biopython parsed");
                continue;
            }

            if (expected.ContainsKey("error"))
            {
                failures.Add($"{r.Id} parsed but Biopython refused ({expected["error"]})");
                continue;
            }

            var actual = new JsonObject
            {
                ["alignments"] = new JsonArray(file.Alignments.Select(BiopythonCoordOracle.ToSpec).ToArray()),
                ["metadata"] = new JsonObject(),
                ["targets"] = new JsonArray(file.Targets.Select(t => BiopythonCoordOracle.ToSpec(new BioCoordAlignment { Sequences = [t], Coordinates = [] })["sequences"]![0]!.DeepClone()).ToArray()),
                ["declaration"] = new JsonObject
                {
                    ["name"] = file.Declaration!.Name,
                    ["comment"] = file.Declaration.Comment,
                    ["fields"] = new JsonArray(file.Declaration.Fields.Select(f => (JsonNode)new JsonArray(f.AsType, f.Name, f.Comment)).ToArray()),
                },
            };
            if (!BiopythonCoordOracle.Same(actual, expected))
                failures.Add($"{r.Id}\n--- expected\n{expected.ToJsonString()[..Math.Min(900, expected.ToJsonString().Length)]}\n--- actual\n{actual.ToJsonString()[..Math.Min(900, actual.ToJsonString().Length)]}");
        }

        Assert.That(failures, Is.Empty, string.Join("\n=====\n", failures.Take(3)) + $"\n({failures.Count})");
    }
}
