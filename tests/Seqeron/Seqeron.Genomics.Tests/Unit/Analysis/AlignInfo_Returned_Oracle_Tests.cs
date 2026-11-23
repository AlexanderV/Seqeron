using System.IO.Compression;
using System.Reflection;
using System.Text.Json;
using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// B09 WP67 (F308, R14-13): <see cref="SummaryInfo"/> / <see cref="Pssm"/> equal Bio.Align.AlignInfo of Biopython 1.79 on every call the user ran (handoff/R14-13-aligninfo-biopython179/returned: 64 alignments, 6 199 calls, Python 3.9.25): the strings of
/// <c>dumb_consensus</c> / <c>gap_consensus</c> / <c>get_column</c> / <c>print_info_content</c>, the list of <c>replacement_dictionary</c>, the axis / rows / <c>str()</c> of <c>pos_specific_score_matrix</c>, the value and the <c>ic_vector</c> of <c>information_content</c> (floats to the
/// bit) and, on every call that raised, the kind and the message of the exception. The one family of calls the typed API cannot express is <c>chars_to_ignore</c> given as a string (Python's <c>isinstance(..., list)</c> check): the test requires the recorded <c>TypeError</c> and counts them.
/// </summary>
[TestFixture]
public class AlignInfo_Returned_Oracle_Tests
{
    private static IEnumerable<JsonElement> Lines(string name)
    {
        string resource = $"Seqeron.Genomics.Tests.TestData.Alignment.{name}";
        using var stream = Assembly.GetExecutingAssembly().GetManifestResourceStream(resource) ?? throw new InvalidOperationException($"Embedded oracle '{resource}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        string? line;
        var all = new List<JsonElement>();
        while ((line = reader.ReadLine()) is not null)
            all.Add(JsonDocument.Parse(line).RootElement.Clone());
        return all;
    }

    private static string Message(Exception e) => e is ArgumentException a && a.ParamName is not null ? a.Message[..a.Message.IndexOf(" (Parameter", StringComparison.Ordinal)] : e.Message;

    private static bool KindMatches(string kind, Exception e) => kind switch
    {
        "ValueError" => e.GetType() == typeof(ArgumentException),
        "IndexError" => e is ArgumentOutOfRangeException or InvalidOperationException,
        "TypeError" or "AssertionError" => e is InvalidOperationException,
        "KeyError" => e is KeyNotFoundException,
        "ZeroDivisionError" => e is DivideByZeroException,
        _ => false,
    };

    private static Dictionary<char, double>? Table(JsonElement kw) => kw.TryGetProperty("e_freq_table", out var t) ? t.EnumerateObject().ToDictionary(p => p.Name[0], p => p.Value.GetDouble()) : null;

    private static List<string>? Ignore(JsonElement kw) => kw.TryGetProperty("chars_to_ignore", out var c) ? [.. c.EnumerateArray().Select(x => x.GetString()!)] : null;

    private static double Num(JsonElement kw, string name, double dflt) => kw.TryGetProperty(name, out var v) ? v.GetDouble() : dflt;

    // the call on a fresh summary: the result in the encoding of run_cases.py
    private static JsonNode? Call(MultipleSeqAlignment alignment, JsonElement call)
    {
        string method = call.GetProperty("method").GetString()!;
        var kw = call.GetProperty("kwargs");
        var summary = new SummaryInfo(alignment);
        switch (method)
        {
            case "dumb_consensus":
                return summary.DumbConsensus(Num(kw, "threshold", 0.7), kw.TryGetProperty("ambiguous", out var a) ? a.GetString()! : "X", kw.TryGetProperty("require_multiple", out var rm) && rm.GetBoolean());
            case "gap_consensus":
                return summary.GapConsensus(Num(kw, "threshold", 0.7), kw.TryGetProperty("ambiguous", out var b) ? b.GetString()! : "X", kw.TryGetProperty("require_multiple", out var rm2) && rm2.GetBoolean());
            case "replacement_dictionary":
            {
                var d = summary.ReplacementDictionary(kw.TryGetProperty("letters", out var l) ? l.GetString() : null, kw.TryGetProperty("skip_chars", out var sc) ? string.Concat(sc.EnumerateArray().Select(x => x.GetString())) : null);
                var arr = new JsonArray();
                foreach (var (key, value) in d)
                    arr.Add((JsonNode)new JsonArray((JsonNode)key.First.ToString(), (JsonNode)key.Second.ToString(), (JsonNode)value));
                return arr;
            }

            case "pos_specific_score_matrix":
            {
                if (kw.TryGetProperty("chars_to_ignore", out var ci) && ci.ValueKind == JsonValueKind.String)
                    throw new NotSupportedException("chars_to_ignore as a string");
                string? axis = kw.TryGetProperty("axis_seq", out var ax) ? (ax.GetString() == "dumb" ? summary.DumbConsensus() : ax.GetString()) : null;
                var p = summary.PosSpecificScoreMatrix(axis, Ignore(kw));
                var rows = new JsonArray();
                var axisArr = new JsonArray();
                for (int i = 0; i < p.Count; i++)
                {
                    axisArr.Add((JsonNode)p.GetResidue(i));
                    var counts = new JsonArray();
                    foreach (var (k, v) in p[i])
                        counts.Add((JsonNode)new JsonArray((JsonNode)k.ToString(), (JsonNode)v));
                    rows.Add((JsonNode)new JsonArray((JsonNode)p.GetResidue(i), counts));
                }

                return new JsonObject { ["axis"] = axisArr, ["rows"] = rows, ["str"] = p.Render() };
            }

            case "information_content":
            case "ic_vector_after":
            {
                double v = summary.InformationContent((int)Num(kw, "start", 0), kw.TryGetProperty("end", out var e) ? (int)e.GetDouble() : null, Table(kw), Num(kw, "log_base", 2), Ignore(kw), Num(kw, "pseudo_count", 0));
                var vec = new JsonArray();
                foreach (double x in summary.IcVector)
                    vec.Add((JsonNode)x);
                return method == "ic_vector_after" ? vec : new JsonObject { ["value"] = v, ["ic_vector"] = vec };
            }

            case "print_info_content":
                if (kw.EnumerateObject().Any())
                    summary.InformationContent(eFreqTable: Table(kw), charsToIgnore: Ignore(kw));
                return summary.PrintInfoContent();
            case "get_column":
                return summary.GetColumn(kw.GetProperty("col").GetInt32());
            default:
                throw new InvalidOperationException(method);
        }
    }

    private static string Describe(JsonNode? n) => n?.ToJsonString() ?? "null";

    private static bool Same(JsonElement expected, JsonNode? got)
    {
        switch (expected.ValueKind)
        {
            case JsonValueKind.String:
                return got is JsonValue sv && sv.TryGetValue<string>(out string? s) && s == expected.GetString();
            case JsonValueKind.Number:
                return got is JsonValue nv && nv.TryGetValue<double>(out double d) && d == expected.GetDouble();
            case JsonValueKind.Array:
            {
                if (got is not JsonArray arr || arr.Count != expected.GetArrayLength())
                    return false;
                int i = 0;
                foreach (var x in expected.EnumerateArray())
                {
                    if (!Same(x, arr[i++]))
                        return false;
                }

                return true;
            }

            case JsonValueKind.Object:
            {
                if (got is not JsonObject obj)
                    return false;
                foreach (var p in expected.EnumerateObject())
                {
                    if (!obj.TryGetPropertyValue(p.Name, out var g) || !Same(p.Value, g))
                        return false;
                }

                return true;
            }

            default:
                return false;
        }
    }

    [Test]
    public void SummaryInfo_EqualsBiopython179_AllReturnedCalls()
    {
        var cases = Lines("r14_13_aligninfo_cases.jsonl.gz").ToList();
        var results = Lines("r14_13_aligninfo_results.jsonl.gz").ToList();
        Assert.That(results, Has.Count.EqualTo(cases.Count));
        var failures = new List<string>();
        int calls = 0, ok = 0, errors = 0, untyped = 0;
        for (int c = 0; c < cases.Count; c++)
        {
            var cs = cases[c];
            string id = cs.GetProperty("id").GetString()!;
            Assert.That(results[c].GetProperty("id").GetString(), Is.EqualTo(id));
            var alignment = new MultipleSeqAlignment(cs.GetProperty("records").EnumerateArray().Select(r =>
            {
                var rec = new SeqRecord(r.GetProperty("seq").GetString()!, r.GetProperty("id").GetString());
                if (r.TryGetProperty("weight", out var w))
                    rec.Annotations["weight"] = w.GetDouble();
                return rec;
            }));
            var callList = cs.GetProperty("calls").EnumerateArray().ToList();
            var expectedList = results[c].GetProperty("results").EnumerateArray().ToList();
            Assert.That(expectedList, Has.Count.EqualTo(callList.Count));
            for (int k = 0; k < callList.Count; k++)
            {
                calls++;
                var call = callList[k];
                var expected = expectedList[k];
                string where = $"{id} call {k} {call.GetProperty("method").GetString()} {call.GetProperty("kwargs")}";
                try
                {
                    var got = Call(alignment, call);
                    if (expected.TryGetProperty("error", out var err))
                        failures.Add($"{where}: no error, expected {err}");
                    else if (!Same(expected.GetProperty("ok"), got))
                        failures.Add($"{where}: got {Describe(got)} expected {expected.GetProperty("ok")}");
                    else
                        ok++;
                }
                catch (NotSupportedException)
                {
                    // chars_to_ignore given as a string: Python checks isinstance(chars_to_ignore, list) at run time; the typed API takes a list
                    if (expected.TryGetProperty("error", out var err) && err[0].GetString() == "TypeError" && err[1].GetString() == "chars_to_ignore should be a list.")
                        untyped++;
                    else
                        failures.Add($"{where}: string chars_to_ignore but expected {expected}");
                }
                catch (Exception e)
                {
                    errors++;
                    if (!expected.TryGetProperty("error", out var err))
                    {
                        failures.Add($"{where}: {e.GetType().Name}: {e.Message} but expected {expected.GetProperty("ok")}");
                        continue;
                    }

                    string kind = err[0].GetString()!, message = err[1].GetString()!;
                    if (!KindMatches(kind, e) || (kind != "AssertionError" && Message(e) != message))
                        failures.Add($"{where}: {e.GetType().Name}: {Message(e)} expected {kind}: {message}");
                }
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {calls}:\n" + string.Join("\n", failures.Take(12)));
        Assert.That(calls, Is.EqualTo(6199));
        Assert.That(ok, Is.GreaterThan(3000), "calls that returned a value, equal to the reference");
        Assert.That(errors, Is.GreaterThan(1500), "calls that raised, with the same kind and message");
        Assert.That(untyped, Is.EqualTo(128), "calls with chars_to_ignore as a string (TypeError of the reference; not expressible with the typed API)");
    }
}
