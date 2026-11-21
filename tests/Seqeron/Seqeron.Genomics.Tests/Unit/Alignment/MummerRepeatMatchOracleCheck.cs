using System.Globalization;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>The oracle comparison of <c>repeat-match</c> and <c>exact-tandems</c> (no test framework: also run from a scratch console).</summary>
public static class MummerRepeatMatchOracleCheck
{
    public static MummerRepeatMatchOptions ParseArgs(IReadOnlyList<string> args)
    {
        var o = new MummerRepeatMatchOptions();
        var rest = new Queue<string>(args);
        while (rest.Count > 0)
        {
            string arg = rest.Dequeue();
            switch (arg)
            {
                case "-f": o = o with { ForwardOnly = true }; break;
                case "-t": o = o with { TandemOnly = true }; break;
                case "-E": o = o with { Exhaustive = true }; break;
                case "-n": o = o with { MinMatchLength = int.Parse(rest.Dequeue(), CultureInfo.InvariantCulture) }; break;
                case "-V": o = o with { Verbose = int.Parse(rest.Dequeue(), CultureInfo.InvariantCulture) }; break;
                default: throw new ArgumentException("unknown option " + arg);
            }
        }

        return o;
    }

    public static (List<string> Failures, Dictionary<string, int> Counts) CheckAll(string resource = "repeatmatch_oracle.jsonl.gz")
    {
        var failures = new List<string>();
        var counts = new Dictionary<string, int>();
        void Count(string key) => counts[key] = counts.GetValueOrDefault(key) + 1;
        int n = 0;
        foreach (var j in MummerOracleData.Rows(resource))
        {
            n++;
            string kind = j.Str("kind");
            int rc = j.GetProperty("rc").GetInt32();
            string label = $"row {n} {kind}";
            if (kind == "rm")
            {
                var args = j.StrList("args");
                label += " " + string.Join(' ', args);
                MummerRepeatMatchResult? r = null;
                string error = "<error>";
                try
                {
                    r = MummerRepeatMatch.Run(j.Str("fasta"), ParseArgs(args));
                }
                catch (Exception e) when (e is FormatException or ArgumentException)
                {
                    error = e.Message;
                }

                if (j.TryGetProperty("ub", out var ub) && ub.GetBoolean())
                {
                    Count("rmUndefined");
                    if (r is not null)
                        failures.Add(label + ": undefined in the program, the port returns text");
                    continue;
                }

                if (rc == 255 && j.Str("out").Length == 0 && j.Str("err").StartsWith("ERROR:  Illegal", StringComparison.Ordinal))
                {
                    Count("rmIllegal");
                    if (r is not null || j.Str("err") != error + "\n")
                        failures.Add(label + ": the illegal option is not refused with the message of the program");
                    continue;
                }

                if (rc != 0 && rc != 255)
                {
                    Count("rmAbort");
                    if (r is not null)
                        failures.Add($"{label}: the program aborts (status {rc}), the port returns text");
                    continue;
                }

                if (r is null)
                    failures.Add($"{label}: {error}");
                else if (r.Text != j.Str("out"))
                    failures.Add($"{label}: stdout {MummerOracleData.Diff(r.Text, j.Str("out"))}");
                else if (r.Messages != j.Str("err"))
                    failures.Add($"{label}: stderr {MummerOracleData.Diff(r.Messages, j.Str("err"))}");
                else if (r.ExitCode != rc)
                    failures.Add($"{label}: exit status {r.ExitCode} expected {rc}");
                else
                    Count(rc == 0 ? "rm" : "rmStops");
            }
            else if (kind is "et" or "etx")
            {
                if (j.TryGetProperty("args", out _))
                {
                    Count("etUsage");
                    continue;
                }

                int minLen = j.GetProperty("n").GetInt32();
                MummerExactTandemsResult? r = null;
                string error = "<error>";
                try
                {
                    r = kind == "et" ? MummerExactTandems.RunInstalled(j.Str("fasta"), minLen) : MummerExactTandems.Run(j.Str("fasta"), minLen);
                }
                catch (Exception e) when (e is FormatException or ArgumentException)
                {
                    error = e.Message;
                }

                if (r is null)
                    failures.Add($"{label} n={minLen}: {error}");
                else if (r.Text != j.Str("out"))
                    failures.Add($"{label} n={minLen}: stdout {MummerOracleData.Diff(r.Text, j.Str("out"))}");
                else if (r.Log != j.Str("err"))
                    failures.Add($"{label} n={minLen}: stderr {MummerOracleData.Diff(r.Log, j.Str("err"))}");
                else if (r.ExitCode != rc)
                    failures.Add($"{label} n={minLen}: exit status {r.ExitCode} expected {rc}");
                else if (kind == "et" && !r.LeftBehind.SequenceEqual(j.StrList("left").Select(f => "$$." + f)))
                    failures.Add($"{label} n={minLen}: files left behind {string.Join(',', r.LeftBehind)} expected {string.Join(',', j.StrList("left"))}");
                else
                    Count(kind);
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
