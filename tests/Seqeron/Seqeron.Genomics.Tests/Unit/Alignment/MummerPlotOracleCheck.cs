using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>The oracle comparison of <c>mummerplot</c> (no test framework: also run from a scratch console).</summary>
public static class MummerPlotOracleCheck
{
    public static MummerPlotOptions ParseArgs(IReadOnlyList<string> args, string matchName)
    {
        var o = new MummerPlotOptions { MatchName = matchName };
        var rest = new Queue<string>(args);
        while (rest.Count > 0)
        {
            string arg = rest.Dequeue();
            switch (arg)
            {
                case "-p": o = o with { Prefix = rest.Dequeue() }; break;
                case "--postscript": o = o with { Postscript = true }; break;
                case "--png": o = o with { Png = true }; break;
                case "--x11": o = o with { X11 = true }; break;
                case "-t": o = o with { Terminal = rest.Dequeue() }; break;
                case "-s": o = o with { Size = rest.Dequeue() }; break;
                case "--small": o = o with { Small = true }; break;
                case "--medium": o = o with { Medium = true }; break;
                case "--large": o = o with { Large = true }; break;
                case "--color": o = o with { Color = true }; break;
                case "--nocolor": o = o with { Color = false }; break;
                case "-c": case "--coverage": o = o with { Coverage = true }; break;
                case "--nocoverage": o = o with { Coverage = false }; break;
                case "-f": o = o with { Filter = true }; break;
                case "-l": o = o with { Layout = true }; break;
                case "--fat": o = o with { Fat = true }; break;
                case "-b": o = o with { BreakLength = int.Parse(rest.Dequeue(), CultureInfo.InvariantCulture) }; break;
                case "-S": o = o with { Snp = true }; break;
                case "-rv": o = o with { ReverseVideo = true }; break;
                case "--title": o = o with { Title = rest.Dequeue() }; break;
                case "-x": o = o with { XRange = rest.Dequeue() }; break;
                case "-y": o = o with { YRange = rest.Dequeue() }; break;
                case "-rport": o = o with { RPort = int.Parse(rest.Dequeue(), CultureInfo.InvariantCulture) }; break;
                case "-qport": o = o with { QPort = int.Parse(rest.Dequeue(), CultureInfo.InvariantCulture) }; break;
                case "-r": o = o with { IdR = rest.Dequeue() }; break;
                case "-q": o = o with { IdQ = rest.Dequeue() }; break;
                case "-R": o = o with { RFile = rest.Dequeue() }; break;
                case "-Q": o = o with { QFile = rest.Dequeue() }; break;
                default: throw new ArgumentException("unknown option " + arg);
            }
        }

        return o;
    }

    /// <summary>The arguments of a command line the way <c>sh</c> splits it (single quotes only).</summary>
    public static List<string> ShellWords(string command)
    {
        var words = new List<string>();
        var cur = new StringBuilder();
        bool inQuote = false, any = false;
        foreach (char c in command)
        {
            if (c == '\'')
            {
                inQuote = !inQuote;
                any = true;
            }
            else if (c == ' ' && !inQuote)
            {
                if (any || cur.Length > 0)
                    words.Add(cur.ToString());
                cur.Clear();
                any = false;
            }
            else
            {
                cur.Append(c);
            }
        }

        if (any || cur.Length > 0)
            words.Add(cur.ToString());
        return words;
    }

    public static (List<string> Failures, Dictionary<string, int> Counts) CheckAll(string resource = "mummerplot_oracle.jsonl.gz")
    {
        var failures = new List<string>();
        var counts = new Dictionary<string, int>();
        void Count(string key) => counts[key] = counts.GetValueOrDefault(key) + 1;
        int n = 0;
        foreach (var j in MummerOracleData.Rows(resource))
        {
            n++;
            var args = j.StrList("args");
            string label = $"row {n} mummerplot {string.Join(' ', args)}";
            if (j.TryGetProperty("ub", out var ub) && ub.GetBoolean())
            {
                Count("mpUndefined");
                continue;
            }

            var files = j.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
            MummerPlotResult r;
            try
            {
                r = MummerPlot.Run(j.Str("match"), files, ParseArgs(args, j.Str("name")));
            }
            catch (Exception e) when (e is FormatException or InvalidOperationException or ArgumentException)
            {
                failures.Add($"{label}: <error: {e.Message}>");
                continue;
            }

            int rc = j.GetProperty("rc").GetInt32();
            var outfiles = j.GetProperty("outfiles").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
            string? cmd = j.GetProperty("cmd").ValueKind == System.Text.Json.JsonValueKind.Null ? null : string.Join('\n', j.GetProperty("cmd").EnumerateArray().Select(x => x.GetString()!));
            string? ours = r.GnuplotCommand is null ? null : string.Join('\n', ShellWords(r.GnuplotCommand).Skip(1));
            if (r.Stdout != j.Str("out"))
                failures.Add($"{label}: stdout {MummerOracleData.Diff(r.Stdout, j.Str("out"))}");
            else if (r.Stderr != j.Str("err"))
                failures.Add($"{label}: stderr {MummerOracleData.Diff(r.Stderr, j.Str("err"))}");
            else if ((rc == 0) != (r.ExitCode == 0))
                failures.Add($"{label}: exit status {r.ExitCode} expected {rc}");
            else if (!outfiles.Keys.OrderBy(k => k, StringComparer.Ordinal).SequenceEqual(r.Files.Keys.OrderBy(k => k, StringComparer.Ordinal)))
                failures.Add($"{label}: files {string.Join(',', r.Files.Keys)} expected {string.Join(',', outfiles.Keys)}");
            else if (ours != cmd)
                failures.Add($"{label}: gnuplot call '{ours}' expected '{cmd}'");
            else
            {
                bool ok = true;
                foreach (var (fn, text) in outfiles)
                {
                    if (r.Files[fn] != text)
                    {
                        failures.Add($"{label}: {fn} {MummerOracleData.Diff(r.Files[fn], text)}");
                        ok = false;
                        break;
                    }
                }

                if (ok)
                    Count(rc == 0 ? "mp" : "mpDies");
            }
        }

        counts["rows"] = n;
        return (failures, counts);
    }
}
