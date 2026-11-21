using System.Globalization;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>The oracle comparison of <c>mapview</c> (no test framework: also run from a scratch console).</summary>
public static class MummerMapviewOracleCheck
{
    public static MummerMapviewOptions ParseArgs(IReadOnlyList<string> args, string coordsName, bool same)
    {
        var o = new MummerMapviewOptions { CoordsName = coordsName, CdsName = same ? "utr.gff" : "cds.gff" };
        var rest = new Queue<string>(args);
        while (rest.Count > 0)
        {
            string arg = rest.Dequeue();
            switch (arg)
            {
                case "-n": o = o with { FileCount = int.Parse(rest.Dequeue(), CultureInfo.InvariantCulture) }; break;
                case "-p": o = o with { Prefix = rest.Dequeue() }; break;
                case "-v": o = o with { Verbose = true }; break;
                case "-d": o = o with { MaxDistance = int.Parse(rest.Dequeue(), CultureInfo.InvariantCulture) }; break;
                case "-I": o = o with { PrintQueryNames = true }; break;
                case "-Ir": o = o with { PrintGeneNames = true }; break;
                case "-f": o = o with { Format = rest.Dequeue() }; break;
                case "-m": o = o with { Magnification = double.Parse(rest.Dequeue(), CultureInfo.InvariantCulture) }; break;
                case "-x1": o = o with { X1 = int.Parse(rest.Dequeue(), CultureInfo.InvariantCulture) }; break;
                case "-x2": o = o with { X2 = int.Parse(rest.Dequeue(), CultureInfo.InvariantCulture) }; break;
                case "-g": o = o with { MgapsReference = rest.Dequeue() }; break;
                default: throw new ArgumentException("unknown option " + arg);
            }
        }

        return o;
    }

    public static (List<string> Failures, Dictionary<string, int> Counts) CheckAll(string resource = "mapview_oracle.jsonl.gz")
    {
        var failures = new List<string>();
        var counts = new Dictionary<string, int>();
        void Count(string key) => counts[key] = counts.GetValueOrDefault(key) + 1;
        int n = 0;
        foreach (var j in MummerOracleData.Rows(resource))
        {
            n++;
            var args = j.StrList("args");
            string label = $"row {n} mapview {string.Join(' ', args)}";
            if (j.TryGetProperty("ub", out var ub) && ub.GetBoolean())
            {
                Count("mvUndefined");
                continue;
            }

            string? utr = j.StrOrNull("utr"), cds = j.StrOrNull("cds");
            bool same = j.GetProperty("same").GetBoolean();
            string name = j.Str("name");
            if (same)
                cds = utr;
            MummerMapviewResult r;
            try
            {
                r = MummerMapview.Run(j.Str("coords"), utr, cds, ParseArgs(args, name, same));
            }
            catch (Exception e) when (e is FormatException or InvalidOperationException or ArgumentException)
            {
                failures.Add($"{label}: <error: {e.Message}>");
                continue;
            }

            var figs = j.GetProperty("figs").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
            var side = j.GetProperty("side").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
            var commands = j.StrList("commands");
            if (r.Stdout != j.Str("out"))
                failures.Add($"{label}: stdout {MummerOracleData.Diff(r.Stdout, j.Str("out"))}");
            else if (r.Stderr != j.Str("err"))
                failures.Add($"{label}: stderr {MummerOracleData.Diff(r.Stderr, j.Str("err"))}");
            else if (r.ExitCode != j.GetProperty("rc").GetInt32())
                failures.Add($"{label}: exit status {r.ExitCode} expected {j.GetProperty("rc").GetInt32()}");
            else if (!figs.Keys.OrderBy(k => k, StringComparer.Ordinal).SequenceEqual(r.Figs.Keys.OrderBy(k => k, StringComparer.Ordinal)))
                failures.Add($"{label}: figs {string.Join(',', r.Figs.Keys)} expected {string.Join(',', figs.Keys)}");
            else if (!r.Commands.SequenceEqual(commands))
                failures.Add($"{label}: commands {string.Join(';', r.Commands)} expected {string.Join(';', commands)}");
            else if (!side.OrderBy(k => k.Key, StringComparer.Ordinal).SequenceEqual(r.Side.OrderBy(k => k.Key, StringComparer.Ordinal)))
                failures.Add($"{label}: side file");
            else
            {
                bool ok = true;
                foreach (var (fn, text) in figs)
                {
                    if (r.Figs[fn] != text)
                    {
                        failures.Add($"{label}: {fn} {MummerOracleData.Diff(r.Figs[fn], text)}");
                        ok = false;
                        break;
                    }
                }

                if (ok)
                    Count(commands.Count > 0 ? "mvConvert" : "mv");
            }
        }

        counts["rows"] = n;
        return (failures, counts);
    }
}
