using System.IO.Compression;
using System.Reflection;
using System.Text.Json;
using System.Text.Json.Nodes;
using System.Text.RegularExpressions;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// B09 WP67 (F307, R14-10b): <see cref="MultipleSeqAlignment"/> / <see cref="SeqRecord"/> equal the legacy Bio.Align.MultipleSeqAlignment of Biopython 1.88 on the rows of TestData/Alignment/multiple_seq_alignment_oracle.jsonl.gz
/// (scripts/alignment/generate_multiple_seq_alignment_oracle.py): every row builds an alignment (or fails as the constructor of the reference) and runs a list of operations (indexing, slicing, add, sort, substitutions, the alignment property,
/// append / extend / del, the column annotations); results and error messages are compared one by one.
/// </summary>
[TestFixture]
public class MultipleSeqAlignment_Oracle_Tests
{
    private static JsonObject Obj(JsonElement e, string name) => e.TryGetProperty(name, out var v) ? (JsonObject)JsonNode.Parse(v.GetRawText())! : [];

    private static SeqRecord Record(JsonElement d) => new(d.GetProperty("seq").GetString()!, d.GetProperty("id").GetString(), d.GetProperty("name").GetString(), d.GetProperty("description").GetString(),
        Obj(d, "annotations"), Obj(d, "letter_annotations"), d.TryGetProperty("dbxrefs", out var x) ? x.EnumerateArray().Select(v => v.GetString()!) : null);

    private static PySlice Slice(JsonElement e)
    {
        int? At(int i) => e[i].ValueKind == JsonValueKind.Null ? null : e[i].GetInt32();
        return new PySlice(At(0), At(1), At(2));
    }

    private static JsonObject Snap(SeqRecord r) => new()
    {
        ["id"] = r.Id, ["name"] = r.Name, ["description"] = r.Description, ["seq"] = r.Seq, ["annotations"] = r.Annotations.DeepClone(), ["letter_annotations"] = r.LetterAnnotations.DeepClone(),
        ["dbxrefs"] = new JsonArray([.. r.DbXrefs.Select(x => (JsonNode)JsonValue.Create(x)!)]),
    };

    private static JsonObject Snap(MultipleSeqAlignment m) => new()
    {
        ["records"] = new JsonArray([.. m.Select(r => (JsonNode)Snap(r))]), ["annotations"] = m.Annotations.DeepClone(), ["column_annotations"] = m.ColumnAnnotations,
    };

    private static string Message(Exception e) => Regex.Replace(e.Message, @" \(Parameter '[^']*'\)$", "");

    private static MultipleSeqAlignment Build(JsonElement row) => new(row.GetProperty("records").EnumerateArray().Select(Record), Obj(row, "annotations"), Obj(row, "column_annotations"));

    private static JsonNode? Run(MultipleSeqAlignment m, JsonElement op)
    {
        switch (op.GetProperty("op").GetString())
        {
            case "len": return new JsonObject { ["kind"] = "int", ["value"] = m.Count };
            case "alignment_length": return new JsonObject { ["kind"] = "int", ["value"] = m.GetAlignmentLength() };
            case "str": return new JsonObject { ["kind"] = "str", ["value"] = m.ToString() };
            case "row": return new JsonObject { ["kind"] = "rec", ["value"] = Snap(m.Row(op.GetProperty("i").GetInt32())) };
            case "rows": return new JsonObject { ["kind"] = "msa", ["value"] = Snap(m.Rows(Slice(op.GetProperty("s")))) };
            case "row_cols": return new JsonObject { ["kind"] = "rec", ["value"] = Snap(m.Row(op.GetProperty("i").GetInt32(), Slice(op.GetProperty("c")))) };
            case "cell": return new JsonObject { ["kind"] = "str", ["value"] = m.Cell(op.GetProperty("i").GetInt32(), op.GetProperty("j").GetInt32()) };
            case "column": return new JsonObject { ["kind"] = "str", ["value"] = m.Column(Slice(op.GetProperty("s")), op.GetProperty("j").GetInt32()) };
            case "slice": return new JsonObject { ["kind"] = "msa", ["value"] = Snap(m.Slice(Slice(op.GetProperty("s")), Slice(op.GetProperty("c")))) };
            case "add":
                var o = op.GetProperty("other");
                return new JsonObject { ["kind"] = "msa", ["value"] = Snap(m.Add(Build(o))) };
            case "sort":
                Func<SeqRecord, IComparable>? key = op.GetProperty("key").GetString() switch
                {
                    "len" => r => r.Length,
                    "seq" => r => r.Seq,
                    "description" => r => r.Description,
                    "gaps" => r => r.Seq.Count(c => c == '-'),
                    _ => null,
                };
                m.Sort(key, op.GetProperty("reverse").GetBoolean());
                return new JsonObject { ["kind"] = "msa", ["value"] = Snap(m) };
            case "substitutions":
                var s = m.Substitutions();
                var values = new JsonArray();
                for (int i = 0; i < s.Letters.Length; i++)
                    values.Add(new JsonArray([.. Enumerable.Range(0, s.Letters.Length).Select(j => (JsonNode)JsonValue.Create(s.Values[i, j])!)]));
                return new JsonObject { ["kind"] = "matrix", ["letters"] = s.Letters, ["values"] = values };
            case "alignment":
                var a = m.ToCoordinateAlignment();
                var rows = new JsonArray();
                foreach (var r in a.Sequences)
                    rows.Add(new JsonObject { ["id"] = r.Id, ["description"] = r.Description, ["seq"] = r.Slice(0, r.Length), ["annotations"] = r.Annotations.DeepClone(), ["letter_annotations"] = r.LetterAnnotations.DeepClone() });
                return new JsonObject
                {
                    ["kind"] = "alignment",
                    ["coordinates"] = new JsonArray([.. a.Coordinates.Select(c => (JsonNode)new JsonArray([.. c.Select(v => (JsonNode)JsonValue.Create(v)!)]))]),
                    ["rows"] = rows,
                    ["annotations"] = a.Attributes["annotations"]!.DeepClone(),
                    ["column_annotations"] = a.Attributes["column_annotations"]!.DeepClone(),
                };
            case "append": m.Append(Record(op.GetProperty("record"))); return new JsonObject { ["kind"] = "msa", ["value"] = Snap(m) };
            case "extend": m.Extend(op.GetProperty("records").EnumerateArray().Select(Record)); return new JsonObject { ["kind"] = "msa", ["value"] = Snap(m) };
            case "del_row": m.RemoveRow(op.GetProperty("i").GetInt32()); return new JsonObject { ["kind"] = "msa", ["value"] = Snap(m) };
            case "del_rows": m.RemoveRows(Slice(op.GetProperty("s"))); return new JsonObject { ["kind"] = "msa", ["value"] = Snap(m) };
            case "set_column_annotation": m.SetColumnAnnotation(op.GetProperty("key").GetString()!, JsonNode.Parse(op.GetProperty("value").GetRawText())); return new JsonObject { ["kind"] = "msa", ["value"] = Snap(m) };
            case "set_column_annotations": m.SetColumnAnnotations((JsonObject)JsonNode.Parse(op.GetProperty("value").GetRawText())!); return new JsonObject { ["kind"] = "msa", ["value"] = Snap(m) };
            case "column_annotations": return new JsonObject { ["kind"] = "msa", ["value"] = Snap(m) };
            default: throw new InvalidOperationException(op.GetProperty("op").GetString());
        }
    }

    // the matrix of substitutions has no kind = alignment-only attributes; the expected value of the reference uses the same names
    private static bool Equal(JsonNode? got, JsonElement expected) => JsonNode.DeepEquals(got, JsonNode.Parse(expected.GetRawText()));

    [Test]
    public void MultipleSeqAlignment_EqualsBiopython_AllRows()
    {
        const string resource = "Seqeron.Genomics.Tests.TestData.Alignment.multiple_seq_alignment_oracle.jsonl.gz";
        using var stream = Assembly.GetExecutingAssembly().GetManifestResourceStream(resource) ?? throw new InvalidOperationException($"Embedded oracle '{resource}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, ops = 0, errors = 0, ctorErrors = 0;
        var kinds = new HashSet<string>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string id = $"row {j.GetProperty("id")}";
            MultipleSeqAlignment m;
            try
            {
                m = Build(j);
            }
            catch (Exception e) when (e is ArgumentException or InvalidOperationException)
            {
                if (!j.TryGetProperty("ctor_error", out var ce))
                    failures.Add($"{id}: constructor failed: {e.Message}");
                else if (Message(e) != ce.GetProperty("message").GetString())
                    failures.Add($"{id}: constructor message '{Message(e)}' expected '{ce.GetProperty("message").GetString()}'");
                ctorErrors++;
                continue;
            }

            if (j.TryGetProperty("ctor_error", out var missing))
            {
                failures.Add($"{id}: constructor succeeded, expected {missing}");
                continue;
            }

            var opList = j.GetProperty("ops").EnumerateArray().ToList();
            var results = j.GetProperty("results").EnumerateArray().ToList();
            for (int k = 0; k < opList.Count; k++)
            {
                ops++;
                var op = opList[k];
                var expected = results[k];
                kinds.Add(op.GetProperty("op").GetString()!);
                string where = $"{id} op {k} {op}";
                JsonNode? got = null;
                Exception? error = null;
                try
                {
                    got = Run(m, op);
                }
                catch (Exception e) when (e is ArgumentException or InvalidOperationException or IndexOutOfRangeException)
                {
                    error = e;
                }

                if (expected.GetProperty("kind").GetString() == "error")
                {
                    errors++;
                    if (error is null)
                        failures.Add($"{where}: no error, expected {expected.GetProperty("message")}");
                    else if (Message(error) != expected.GetProperty("message").GetString())
                        failures.Add($"{where}: message '{Message(error)}' expected '{expected.GetProperty("message").GetString()}'");
                }
                else if (error is not null)
                {
                    failures.Add($"{where}: {error.GetType().Name} {error.Message}");
                }
                else if (!Equal(got, expected))
                {
                    failures.Add($"{where}: got {got?.ToJsonString()} expected {expected}");
                }
            }
        }

        Assert.That(failures, Is.Empty, string.Join(Environment.NewLine, failures.Take(10)));
        Assert.Multiple(() =>
        {
            Assert.That(rows, Is.GreaterThanOrEqualTo(1200));
            Assert.That(ops, Is.GreaterThan(6000));
            Assert.That(errors, Is.GreaterThan(500));
            Assert.That(ctorErrors, Is.GreaterThan(100));
            Assert.That(kinds, Has.Count.GreaterThanOrEqualTo(18), "every operation is exercised");
        });
    }

    [Test]
    public void Add_And_Slice_FromBiopythonDocumentation()
    {
        // the doctests of MultipleSeqAlignment.__add__ / __getitem__ / sort / substitutions
        var left = new MultipleSeqAlignment(
            [new SeqRecord("AAAAC", "Alpha"), new SeqRecord("AAA-C", "Beta"), new SeqRecord("AAAAG", "Gamma")],
            new JsonObject { ["tool"] = "demo", ["name"] = "start" }, new JsonObject { ["stats"] = "CCCXC" });
        var right = new MultipleSeqAlignment(
            [new SeqRecord("GT", "Alpha"), new SeqRecord("GT", "Beta"), new SeqRecord("GT", "Gamma")],
            new JsonObject { ["tool"] = "demo", ["name"] = "end" }, new JsonObject { ["stats"] = "CC" });
        var combined = left.Add(right);
        Assert.Multiple(() =>
        {
            Assert.That(combined.ToString(), Is.EqualTo("Alignment with 3 rows and 7 columns\nAAAACGT Alpha\nAAA-CGT Beta\nAAAAGGT Gamma"));
            Assert.That(combined.Annotations.ToJsonString(), Is.EqualTo("{\"tool\":\"demo\"}"));
            Assert.That(combined.ColumnAnnotations.ToJsonString(), Is.EqualTo("{\"stats\":\"CCCXCCC\"}"));
        });

        var align = new MultipleSeqAlignment([new SeqRecord("AAAACGT", "Alpha"), new SeqRecord("AAA-CGT", "Beta"), new SeqRecord("AAAAGGT", "Gamma"), new SeqRecord("AAAACGT", "Delta"), new SeqRecord("AAA-GGT", "Epsilon")]);
        Assert.Multiple(() =>
        {
            Assert.That(align.Cell(3, 4), Is.EqualTo("C"));
            Assert.That(align.Column(4), Is.EqualTo("CCGCG"));
            Assert.That(align.Column(new PySlice(1, 3), 4), Is.EqualTo("CG"));
            Assert.That(align.Slice(new PySlice(1, 5), new PySlice(3, 6)).ToString(), Is.EqualTo("Alignment with 4 rows and 3 columns\n-CG Beta\nAGG Gamma\nACG Delta\n-GG Epsilon"));
            Assert.That(align.Rows(new PySlice(Step: -1)).Select(r => r.Id), Is.EqualTo(new[] { "Epsilon", "Delta", "Gamma", "Beta", "Alpha" }));
        });

        var substitutions = new MultipleSeqAlignment([new SeqRecord("ACGT", "seq1"), new SeqRecord("A--A", "seq2"), new SeqRecord("ACGT", "seq3"), new SeqRecord("TTTC", "seq4")]).Substitutions();
        Assert.Multiple(() =>
        {
            Assert.That(substitutions.Letters, Is.EqualTo("ACGT"));
            Assert.That(substitutions['A', 'A'], Is.EqualTo(3.0));
            Assert.That(substitutions['A', 'C'], Is.EqualTo(0.5));
            Assert.That(substitutions['A', 'T'], Is.EqualTo(2.5));
            Assert.That(substitutions['C', 'T'], Is.EqualTo(2.0));
            Assert.That(substitutions['G', 'T'], Is.EqualTo(1.0));
        });
    }
}
