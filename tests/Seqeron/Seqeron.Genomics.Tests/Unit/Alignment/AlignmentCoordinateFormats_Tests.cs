using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-27c, B09 F224 / F226): the coordinate formats of Biopython 1.88 <c>Bio.Align</c> (<c>bed</c>, <c>psl</c>, <c>chain</c>, <c>maf</c>): <see cref="AlignmentCoordinateFormats.Write"/> byte for byte and
/// <see cref="AlignmentCoordinateFormats.Parse"/> equal to <c>Bio.Align.parse</c> (scripts/alignment/generate_biopython_coord_oracle.py, biopython_coord_oracle.tsv).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class AlignmentCoordinateFormats_Tests
{
    private static readonly Lazy<List<BiopythonCoordOracle.Row>> Oracle = new(BiopythonCoordOracle.Load);

    private static readonly string[] Formats = ["bed", "psl", "chain", "maf", "mauve", "sam", "exonerate"];

    private static CoordinateWriteOptions Options(JsonObject o) => new()
    {
        BedN = o["bedN"]?.GetValue<int>() ?? 12,
        Header = o["header"]?.GetValue<bool>() ?? true,
        Mask = o["mask"]?.GetValue<string>(),
        Wildcard = o.ContainsKey("wildcard") ? o["wildcard"]?.GetValue<string>() : "N",
        Metadata = o["metadata"] as JsonObject,
        Identifiers = o["identifiers"] is JsonArray ids ? [.. ids.Select(x => x!.GetValue<string>())] : null,
        Md = o["md"]?.GetValue<bool>() ?? false,
        Targets = o["targets"] is JsonArray tg ? [.. tg.Select(t => AlignmentCoordinateFormats.FromJson(new JsonObject { ["sequences"] = new JsonArray(t!.DeepClone()), ["coordinates"] = new JsonArray(), ["attrs"] = new JsonObject() }).Sequences[0])] : null,
        ExonerateFormat = o["exfmt"]?.GetValue<string>() ?? "vulgar",
    };

    [Test]
    public void Oracle_Has200WriterAndReaderRowsPerFormat()
    {
        Assert.Multiple(() =>
        {
            foreach (var f in Formats)
            {
                Assert.That(Oracle.Value.Count(r => r.Format == f && r.Kind == "write"), Is.GreaterThanOrEqualTo(100), f + " writers");
                Assert.That(Oracle.Value.Count(r => r.Format == f && r.Kind == "read"), Is.GreaterThanOrEqualTo(70), f + " readers");
                Assert.That(Oracle.Value.Count(r => r.Format == f), Is.GreaterThanOrEqualTo(200), f);
            }
        });
    }

    [Test]
    public void Write_EveryRow_IsByteEqualToBiopython()
    {
        var failures = new List<string>();
        foreach (var r in Oracle.Value.Where(r => r.Kind == "write" && Formats.Contains(r.Format)))
        {
            AlignmentCoordinateFormats.TryParseFormat(r.Format, out var fmt);
            var aligns = JsonNode.Parse(r.Data1)!.AsArray().Select(s => BiopythonCoordOracle.FromSpec(s!)).ToList();
            string actual;
            try
            {
                actual = AlignmentCoordinateFormats.Write(aligns, fmt, Options(r.Options));
            }
            catch (ArgumentException)
            {
                if (!r.Data2.StartsWith("ERROR:", StringComparison.Ordinal))
                    failures.Add($"{r.Id} {r.Format} threw but Biopython wrote: {r.Data2}");
                continue;
            }

            if (r.Data2.StartsWith("ERROR:", StringComparison.Ordinal))
                failures.Add($"{r.Id} {r.Format} wrote but Biopython refused: {r.Data2}");
            else if (actual != r.Data2)
                failures.Add($"{r.Id} {r.Format} {r.Options}\n--- expected\n{r.Data2}--- actual\n{actual}");
        }

        Assert.That(failures, Is.Empty, string.Join("\n=====\n", failures.Take(4)) + $"\n({failures.Count})");
    }

    [Test]
    public void Parse_EveryRow_EqualsBioAlignParse()
    {
        var failures = new List<string>();
        foreach (var r in Oracle.Value.Where(r => r.Kind == "read" && Formats.Contains(r.Format)))
        {
            AlignmentCoordinateFormats.TryParseFormat(r.Format, out var fmt);
            var expected = JsonNode.Parse(r.Data2)!;
            CoordinateAlignmentFile actualFile;
            try
            {
                actualFile = AlignmentCoordinateFormats.ParseFile(r.Data1, fmt);
            }
            catch (ArgumentException)
            {
                if (!(expected is JsonObject refusedNode && refusedNode.ContainsKey("error")))
                    failures.Add($"{r.Id} {r.Format} threw but Biopython parsed: {r.Data2}");
                continue;
            }

            if (expected is JsonObject err && err.ContainsKey("error"))
            {
                failures.Add($"{r.Id} {r.Format} parsed but Biopython refused ({err["error"]})\n{r.Data1}");
                continue;
            }

            var actualSpec = new JsonArray(actualFile.Alignments.Select(BiopythonCoordOracle.ToSpec).ToArray());
            JsonNode actualNode = actualSpec;
            if (expected is JsonObject fileLevel && fileLevel.ContainsKey("alignments"))
            {
                actualNode = new JsonObject
                {
                    ["alignments"] = actualSpec,
                    ["metadata"] = actualFile.Metadata.DeepClone(),
                    ["targets"] = new JsonArray(actualFile.Targets.Select(t => BiopythonCoordOracle.ToSpec(new BioCoordAlignment { Sequences = [t], Coordinates = [] })["sequences"]![0]!.DeepClone()).ToArray()),
                };
            }

            if (!BiopythonCoordOracle.Same(actualNode, expected))
                failures.Add($"{r.Id} {r.Format}\n--- text\n{r.Data1}--- expected\n{expected.ToJsonString()}\n--- actual\n{actualNode.ToJsonString()}");
        }

        Assert.That(failures, Is.Empty, string.Join("\n=====\n", failures.Take(4)) + $"\n({failures.Count})");
    }
}
