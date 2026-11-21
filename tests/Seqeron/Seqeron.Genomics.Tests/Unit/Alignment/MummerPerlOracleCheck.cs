using System.Globalization;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>The oracle comparison of the perl programs of MUMmer 3.23 (no test framework: also run from a scratch console).</summary>
public static class MummerPerlOracleCheck
{
    public static (List<string> Failures, Dictionary<string, int> Counts) CheckNucmer2Xfig(string resource = "nucmer2xfig_oracle.jsonl.gz")
    {
        var failures = new List<string>();
        int rows = 0, equal = 0, labels = 0, ticks = 0, malformed = 0;
        foreach (var j in MummerOracleData.Rows(resource))
        {
            rows++;
            string got = MummerNucmer2Xfig.Run(j.Str("coords"));
            if (got != j.Str("out"))
            {
                failures.Add($"row {rows}: {MummerOracleData.Diff(got, j.Str("out"))}");
                continue;
            }

            equal++;
            if (got.Contains("4.7124", StringComparison.Ordinal)) labels++;
            if (got.Contains("2 1 0 1 0 7 50", StringComparison.Ordinal)) ticks++;
            if (!j.Str("coords").Contains('|', StringComparison.Ordinal)) malformed++;
        }

        return (failures, new Dictionary<string, int> { ["rows"] = rows, ["equal"] = equal, ["labels"] = labels, ["ticks"] = ticks, ["malformed"] = malformed });
    }

    public const string PrefixWorkDir = "/tmp/b09pfx";

    /// <summary><c>nucmer -p|--prefix</c>: the files and the standard error of the script on every option family.</summary>
    public static (List<string> Failures, Dictionary<string, int> Counts) CheckNucmerPrefix(string resource = "nucmer_prefix_oracle.jsonl.gz")
    {
        var failures = new List<string>();
        int rows = 0, equal = 0, failing = 0, custom = 0, nodelta = 0, coords = 0, subdir = 0;
        foreach (var j in MummerOracleData.Rows(resource))
        {
            rows++;
            var args = j.StrList("args");
            string label = $"row {rows} [{string.Join(' ', args)}]";
            var o = new NucmerOptions();
            using (var e = args.GetEnumerator())
            {
                while (e.MoveNext())
                {
                    switch (e.Current)
                    {
                        case "-p" or "--prefix": e.MoveNext(); o = o with { Prefix = e.Current }; break;
                        case "-l": e.MoveNext(); o = o with { MinMatch = int.Parse(e.Current, CultureInfo.InvariantCulture) }; break;
                        case "-c": e.MoveNext(); o = o with { MinCluster = long.Parse(e.Current, CultureInfo.InvariantCulture) }; break;
                        case "-b": e.MoveNext(); o = o with { BreakLength = int.Parse(e.Current, CultureInfo.InvariantCulture) }; break;
                        case "--nodelta": o = o with { Delta = false }; break;
                        case "-o" or "--coords": o = o with { Coords = true }; break;
                        case "--mum": o = o with { MatchMode = CommonRegionMode.MaximalUnique }; break;
                        case "--maxmatch": o = o with { MatchMode = CommonRegionMode.MaximalExactMatches }; break;
                        case "--mumreference": o = o with { MatchMode = CommonRegionMode.MaximalUniqueReference }; break;
                        case "-f": o = o with { ForwardOnly = true }; break;
                        default: throw new ArgumentException(e.Current);
                    }
                }
            }

            NucmerResult? r = null;
            try
            {
                r = GenomicAnalyzer.Nucmer(j.Str("ref"), j.Str("qry"), PrefixWorkDir + "/ref.fa", PrefixWorkDir + "/qry.fa", o);
            }
            catch (FormatException)
            {
            }

            if (j.GetProperty("rc").GetInt32() != 0)
            {
                // a prefix in a directory that does not exist or two matching algorithms are decided by the script's file and option checks, which an in-memory API has not got
                bool scriptOnly = args.Any(a => a.StartsWith("nodir/", StringComparison.Ordinal)) || (args.Contains("--mum") && args.Contains("--maxmatch"));
                if (!scriptOnly)
                {
                    failing++;
                    if (r is not null)
                        failures.Add(label + ": the script fails, the port returns files");
                }

                continue;
            }

            if (r is null)
            {
                failures.Add(label + ": <error>");
                continue;
            }

            var expected = j.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
            if (!r.Files.Keys.OrderBy(k => k, StringComparer.Ordinal).SequenceEqual(expected.Keys.OrderBy(k => k, StringComparer.Ordinal), StringComparer.Ordinal))
            {
                failures.Add($"{label}: files {string.Join(',', r.Files.Keys)} expected {string.Join(',', expected.Keys)}");
                continue;
            }

            bool ok = true;
            foreach (var (name, text) in r.Files)
            {
                if (text != expected[name])
                {
                    failures.Add($"{label}: {name} {MummerOracleData.Diff(text, expected[name])}");
                    ok = false;
                }
            }

            if (r.Log.TrimEnd('\n') != j.Str("err").Replace("/W/", PrefixWorkDir + "/", StringComparison.Ordinal).TrimEnd('\n'))
            {
                failures.Add($"{label}: stderr {MummerOracleData.Diff(r.Log, j.Str("err").Replace("/W/", PrefixWorkDir + "/", StringComparison.Ordinal))}");
                ok = false;
            }

            if (!ok)
                continue;
            equal++;
            if (o.Prefix != "out") custom++;
            if (!o.Delta) nodelta++;
            if (o.Coords) coords++;
            if (o.Prefix.Contains('/', StringComparison.Ordinal)) subdir++;
        }

        return (failures, new Dictionary<string, int> { ["rows"] = rows, ["equal"] = equal, ["failing"] = failing, ["customPrefix"] = custom, ["nodelta"] = nodelta, ["coords"] = coords, ["subdir"] = subdir });
    }
}
