using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-27b, B09 F126): the Biopython 1.88 <c>Bio.Align</c> writers and readers of <see cref="AlignmentFormats"/> (fasta, a2m, clustal, phylip, nexus, stockholm) against the installed Biopython:
/// TestData/Alignment/biopython_msa_oracle.tsv (scripts/alignment/generate_biopython_msa_oracle.py): every file written byte for byte as <c>Bio.Align.write</c>, and the files Biopython itself can read parsed to the same alignments
/// (ids, descriptions, printed rows, annotations) as <c>Bio.Align.parse</c>.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class AlignmentFormats_Tests
{
    private sealed record Row(string Id, string Format, JsonNode Options, JsonArray Alignments, string File, JsonNode? Parsed);

    private static string Unescape(string s)
    {
        var sb = new System.Text.StringBuilder();
        bool escaped = false;
        foreach (char c in s)
        {
            if (escaped)
            {
                sb.Append(c switch { 'n' => '\n', 't' => '\t', _ => c });
                escaped = false;
            }
            else if (c == '\\')
            {
                escaped = true;
            }
            else
            {
                sb.Append(c);
            }
        }

        return sb.ToString();
    }

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.biopython_msa_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            string parsed = Unescape(f[5]);
            list.Add(new Row(f[0], f[1], JsonNode.Parse(f[2])!, JsonNode.Parse(f[3])!.AsArray(), Unescape(f[4]), parsed == "null" ? null : JsonNode.Parse(parsed)));
        }

        return list;
    });

    private static BioAlignment FromJson(JsonNode node)
    {
        var o = node.AsObject();
        var records = o["records"]!.AsArray().Select(r =>
        {
            var ro = r!.AsObject();
            return new BioAlignmentRecord
            {
                Id = ro["id"]!.GetValue<string>(),
                Description = ro["desc"]!.GetValue<string>(),
                Row = ro["row"]!.GetValue<string>(),
                Annotations = (JsonObject)ro["ann"]!.DeepClone(),
                LetterAnnotations = (JsonObject)ro["letter"]!.DeepClone(),
                DbXrefs = [.. ro["dbx"]!.AsArray().Select(x => x!.GetValue<string>())],
            };
        }).ToList();
        return new BioAlignment { Records = records, Annotations = (JsonObject)o["ann"]!.DeepClone(), ColumnAnnotations = (JsonObject)o["col"]!.DeepClone() };
    }

    private static JsonNode ToJson(BioAlignment a) => new JsonObject
    {
        ["records"] = new JsonArray([.. a.Records.Select(r => (JsonNode)new JsonObject
        {
            ["id"] = r.Id,
            ["desc"] = r.Description,
            ["row"] = r.Row,
            ["ann"] = r.Annotations.DeepClone(),
            ["letter"] = r.LetterAnnotations.DeepClone(),
            ["dbx"] = new JsonArray([.. r.DbXrefs.Select(x => (JsonNode)JsonValue.Create(x)!)]),
        })]),
        ["ann"] = a.Annotations.DeepClone(),
        ["col"] = a.ColumnAnnotations.DeepClone(),
    };

    private static string FirstDifference(string got, string expected)
    {
        var g = got.Split('\n');
        var e = expected.Split('\n');
        for (int i = 0; i < Math.Max(g.Length, e.Length); i++)
        {
            string x = i < g.Length ? g[i] : "<end>";
            string y = i < e.Length ? e[i] : "<end>";
            if (x != y)
                return $"line {i + 1}: got '{x}' expected '{y}'";
        }

        return "";
    }

    private static void CheckWriters(string format)
    {
        Assert.That(AlignmentFormats.TryParseFormat(format, out var fmt), format);
        var failures = new List<string>();
        var rows = Oracle.Value.Where(r => r.Format == format).ToList();
        foreach (var r in rows)
        {
            try
            {
                bool? interleave = r.Options is JsonObject oo && oo.TryGetPropertyValue("interleave", out var iv) && iv is not null ? iv.GetValue<bool>() : null;
                string got = AlignmentFormats.Write([.. r.Alignments.Select(a => FromJson(a!))], fmt, new AlignmentWriteOptions { NexusInterleave = interleave });
                if (got != r.File)
                    failures.Add($"#{r.Id}: {FirstDifference(got, r.File)}");
            }
            catch (Exception e)
            {
                failures.Add($"#{r.Id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(100));
        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}\n" + string.Join("\n", failures.Take(15)));
    }

    private static string Clip(string s) => s.Length <= 300 ? s : s[..300];

    private static void CheckReaders(string format)
    {
        Assert.That(AlignmentFormats.TryParseFormat(format, out var fmt), format);
        var failures = new List<string>();
        var rows = Oracle.Value.Where(r => r.Format == format && r.Parsed is not null).ToList();
        foreach (var r in rows)
        {
            try
            {
                var got = new JsonArray([.. AlignmentFormats.Parse(r.File, fmt).Select(ToJson)]);
                if (!JsonNode.DeepEquals(got, r.Parsed))
                    failures.Add($"#{r.Id}: got {Clip(got.ToJsonString())}\n   expected {Clip(r.Parsed!.ToJsonString())}");
            }
            catch (Exception e)
            {
                failures.Add($"#{r.Id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(100));
        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}\n" + string.Join("\n", failures.Take(10)));
    }

    [Test]
    public void FastaWriter_IsByteIdenticalToBiopython() => CheckWriters("fasta");

    [Test]
    public void A2mWriter_IsByteIdenticalToBiopython() => CheckWriters("a2m");

    [Test]
    public void ClustalWriter_IsByteIdenticalToBiopython() => CheckWriters("clustal");

    [Test]
    public void PhylipWriter_IsByteIdenticalToBiopython() => CheckWriters("phylip");

    [Test]
    public void NexusWriter_IsByteIdenticalToBiopython() => CheckWriters("nexus");

    [Test]
    public void StockholmWriter_IsByteIdenticalToBiopython() => CheckWriters("stockholm");

    [Test]
    public void FastaReader_ParsesLikeBiopython() => CheckReaders("fasta");

    [Test]
    public void A2mReader_ParsesLikeBiopython() => CheckReaders("a2m");

    [Test]
    public void ClustalReader_ParsesLikeBiopython() => CheckReaders("clustal");

    [Test]
    public void PhylipReader_ParsesLikeBiopython() => CheckReaders("phylip");

    [Test]
    public void NexusReader_ParsesLikeBiopython() => CheckReaders("nexus");

    [Test]
    public void StockholmReader_ParsesLikeBiopython() => CheckReaders("stockholm");

    [Test]
    public void UnknownFormatName_IsRejected()
    {
        Assert.That(AlignmentFormats.TryParseFormat("maf", out _), Is.False);
        Assert.That(AlignmentFormats.TryParseFormat(" FASTA ", out var f), Is.True);
        Assert.That(f, Is.EqualTo(AlignmentFileFormat.Fasta));
    }
}
