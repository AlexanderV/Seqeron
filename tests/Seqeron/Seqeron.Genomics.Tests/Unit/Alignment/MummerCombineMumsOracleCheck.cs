using System.Globalization;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>The oracle comparison of <c>combineMUMs</c> and the script <c>run-mummer3</c> (no test framework: also run from a scratch console).</summary>
public static class MummerCombineMumsOracleCheck
{
    public static MummerCombineMumsOptions ParseArgs(IReadOnlyList<string> args, string referenceName = "ref.fa", string queryName = "qry.fa")
    {
        var o = new MummerCombineMumsOptions { ReferenceName = referenceName, QueryName = queryName };
        var rest = new Queue<string>(args);
        while (rest.Count > 0)
        {
            string arg = rest.Dequeue();
            switch (arg)
            {
                case "-D": o = o with { OnlyDifferencePositions = true }; break;
                case "-S": o = o with { ShowDifferences = true }; break;
                case "-n": o = o with { NucleotidesOnly = true }; break;
                case "-x": o = o with { NoCoverFiles = true }; break;
                case "-t": o = o with { TagFromFastaLine = true }; break;
                case "-e": o = o with { ErrorRate = double.Parse(rest.Dequeue(), CultureInfo.InvariantCulture) }; break;
                case "-N": o = o with { ConsecutiveNonAcgt = int.Parse(rest.Dequeue(), CultureInfo.InvariantCulture) }; break;
                case "-v": o = o with { Verbose = int.Parse(rest.Dequeue(), CultureInfo.InvariantCulture) }; break;
                case "-q": o = o with { QuerySuffix = rest.Dequeue() }; break;
                case "-r": o = o with { ReferenceSuffix = rest.Dequeue() }; break;
                case "-W": o = o with { ErrorsFileName = rest.Dequeue() }; break;
                default: throw new ArgumentException("unknown option " + arg);
            }
        }

        return o;
    }

    public static (List<string> Failures, Dictionary<string, int> Counts) CheckAll(string resource = "combinemums_oracle.jsonl.gz")
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
            bool ub = j.TryGetProperty("ub", out var u) && u.ValueKind == System.Text.Json.JsonValueKind.True;
            int rc = j.GetProperty("rc").GetInt32();
            if (kind == "cmums")
            {
                var args = j.StrList("args");
                MummerCombineMumsResult? r = null;
                string error = "<error>";
                try
                {
                    r = MummerCombineMums.Run(j.Str("ref"), j.Str("qry"), j.Str("gaps"), ParseArgs(args));
                }
                catch (FormatException e)
                {
                    error = "<error: " + e.Message + ">";
                }

                if (ub)
                {
                    Count("cmumsUndefined");
                    continue;
                }

                if (rc != 0)
                {
                    Count("cmumsFail");
                    if (r is not null)
                        failures.Add($"{label} {string.Join(' ', args)}: the program fails (status {rc}), the port returns text");
                    continue;
                }

                if (r is null)
                    failures.Add($"{label} {string.Join(' ', args)}: {error}");
                else if (r.Text != j.Str("out"))
                    failures.Add($"{label} {string.Join(' ', args)}: stdout {MummerOracleData.Diff(r.Text, j.Str("out"))}");
                else if (r.Log != j.Str("err"))
                    failures.Add($"{label} {string.Join(' ', args)}: stderr {MummerOracleData.Diff(r.Log, j.Str("err"))}");
                else
                {
                    var files = j.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
                    if (!files.Keys.OrderBy(k => k, StringComparer.Ordinal).SequenceEqual(r.Files.Keys.OrderBy(k => k, StringComparer.Ordinal)))
                        failures.Add($"{label} {string.Join(' ', args)}: files {string.Join(',', r.Files.Keys)} expected {string.Join(',', files.Keys)}");
                    else
                    {
                        bool ok = true;
                        foreach (var (name, text) in files)
                        {
                            if (r.Files[name] != text)
                            {
                                failures.Add($"{label} {string.Join(' ', args)}: {name} {MummerOracleData.Diff(r.Files[name], text)}");
                                ok = false;
                            }
                        }

                        if (ok)
                            Count("cmums");
                    }
                }
            }
            else if (kind == "runmummer3")
            {
                RunMummer3Result? r = null;
                try
                {
                    r = MummerRunPrograms.RunMummer3(j.Str("ref"), j.Str("qry"), "ref.fa", "qry.fa", j.Str("prefix"));
                }
                catch (FormatException)
                {
                }

                if (ub)
                {
                    Count("runmummer3Undefined");
                    continue;
                }

                if (rc != 0)
                {
                    Count("runmummer3Fail");
                    if (r is not null)
                        failures.Add(label + ": the script fails, the port returns files");
                    continue;
                }

                if (r is null)
                {
                    failures.Add(label + ": <error>");
                    continue;
                }

                bool ok = true;
                foreach (var (ext, text) in new[] { ("out", r.Out), ("gaps", r.Gaps), ("align", r.Align), ("errorsgaps", r.ErrorsGaps) })
                {
                    if (text != j.GetProperty("files").GetProperty(ext).GetString())
                    {
                        failures.Add($"{label}: .{ext} {MummerOracleData.Diff(text, j.GetProperty("files").GetProperty(ext).GetString() ?? "<null>")}");
                        ok = false;
                    }
                }

                string err = string.Join('\n', r.Log.Split('\n').Where(l => !l.StartsWith("# CONSTRUCTIONTIME", StringComparison.Ordinal)));
                if (r.Stdout != j.Str("out"))
                {
                    failures.Add($"{label}: stdout {MummerOracleData.Diff(r.Stdout, j.Str("out"))}");
                    ok = false;
                }

                if (err != j.Str("err"))
                {
                    failures.Add($"{label}: stderr {MummerOracleData.Diff(err, j.Str("err"))}");
                    ok = false;
                }

                if (ok)
                    Count("runmummer3");
            }
            else
            {
                failures.Add(label + ": unknown row");
            }
        }

        counts["rows"] = n;
        return (failures, counts);
    }
}
