using System.Text.Json;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>The oracle comparison of <c>gaps</c>, <c>delta2blocks</c>, <c>mummer-annotate</c> and <c>run-mummer1</c> (no test framework: also run from a scratch console).</summary>
public static class MummerGapsOracleCheck
{
    public static (List<string> Failures, Dictionary<string, int> Counts) CheckAll(string resource = "mummer_gaps_oracle.jsonl.gz")
    {
        var failures = new List<string>();
        var counts = new Dictionary<string, int>();
        void Count(string key) => counts[key] = counts.GetValueOrDefault(key) + 1;
        int n = 0;
        foreach (var j in MummerOracleData.Rows(resource))
        {
            n++;
            string kind = j.Str("kind");
            string label = $"row {n} {kind}";
            switch (kind)
            {
                case "gaps":
                {
                    string? got = null;
                    try
                    {
                        got = MummerGaps.Run(j.Str("file"), j.GetProperty("rflag").GetBoolean(), j.Str("input"));
                    }
                    catch (FormatException)
                    {
                    }

                    if (j.GetProperty("rc").GetInt32() != 0)
                    {
                        Count("gapsCrash");
                        if (got is not null)
                            failures.Add(label + ": the program crashes, the port returns text");
                    }
                    else if (got != j.Str("out"))
                        failures.Add($"{label}: {MummerOracleData.Diff(got ?? "<error>", j.Str("out"))}");
                    else
                        Count("gaps");
                    break;
                }

                case "delta2blocks":
                {
                    var args = j.StrList("args");
                    var o = new MummerShowAlignsOptions { SortByQuery = args.Contains("-q"), SortByReference = args.Contains("-r") };
                    string? got = null;
                    bool crashed = false;
                    try
                    {
                        got = MummerDelta2Blocks.Run(j.Str("delta"), j.Str("ref"), j.Str("qry"), args[^2], args[^1], o).Text;
                    }
                    catch (InvalidOperationException)
                    {
                        crashed = true;
                    }

                    if (j.GetProperty("rc").GetInt32() != 0)
                    {
                        Count("delta2blocksCrash");
                        if (!crashed)
                            failures.Add(label + ": the program crashes, the port does not");
                    }
                    else if (got != j.Str("out"))
                        failures.Add($"{label} {string.Join(' ', args)}: {MummerOracleData.Diff(got ?? "<error>", j.Str("out"))}");
                    else
                        Count("delta2blocks");
                    break;
                }

                case "annotate":
                {
                    if (j.GetProperty("ub").GetBoolean())
                    {
                        Count("annotateUndefined");
                        break;
                    }

                    var files = j.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
                    MummerAnnotateResult? r = null;
                    try
                    {
                        r = MummerAnnotate.Run(j.Str("gaps"), j.Str("data"), files);
                    }
                    catch (InvalidOperationException)
                    {
                    }

                    if (r is null)
                        failures.Add(label + ": <error>");
                    else if (r.Text != j.Str("out"))
                        failures.Add($"{label}: stdout {MummerOracleData.Diff(r.Text, j.Str("out"))}");
                    else if (r.WithErrors != j.Str("errors"))
                        failures.Add($"{label}: witherrors.gaps {MummerOracleData.Diff(r.WithErrors, j.Str("errors"))}");
                    else
                        Count("annotate");
                    break;
                }

                case "runmummer1":
                {
                    var args = j.StrList("args");
                    RunMummer1Result? r = null;
                    try
                    {
                        r = MummerRunPrograms.RunMummer1(j.Str("ref"), j.Str("qry"), "ref.fa", j.Str("prefix"), args.Contains("-r"));
                    }
                    catch (Exception e) when (e is FormatException or InvalidOperationException)
                    {
                    }

                    if (j.GetProperty("ub").GetBoolean())
                    {
                        Count("runmummer1Undefined");
                        break;
                    }

                    if (j.GetProperty("rc").GetInt32() != 0)
                    {
                        Count("runmummer1Fail");
                        if (r is not null)
                            failures.Add(label + ": the script fails, the port returns files");
                        break;
                    }

                    if (r is null)
                    {
                        failures.Add(label + ": <error>");
                        break;
                    }

                    bool ok = true;
                    foreach (var (ext, text) in new[] { ("out", r.Out), ("gaps", r.Gaps), ("align", r.Align), ("errorsgaps", r.ErrorsGaps) })
                    {
                        if (text != j.GetProperty("files").GetProperty(ext).GetString())
                        {
                            failures.Add($"{label} {string.Join(' ', args)}: .{ext} {MummerOracleData.Diff(text, j.GetProperty("files").GetProperty(ext).GetString() ?? "<null>")}");
                            ok = false;
                        }
                    }

                    if (ok)
                        Count("runmummer1");
                    break;
                }

                default:
                    failures.Add(label + ": unknown row");
                    break;
            }
        }

        counts["rows"] = n;
        return (failures, counts);
    }
}
