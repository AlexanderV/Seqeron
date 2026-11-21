using System.Globalization;
using System.Text.Json;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>The oracle comparison of <c>promer</c> and its stages (no test framework: also run from a scratch console).</summary>
public static class MummerPromerOracleCheck
{
    public const string WorkDir = "/tmp/b09pro";

    internal static PromerOptions ParseArgs(IReadOnlyList<string> args)
    {
        var o = new PromerOptions();
        using var e = args.GetEnumerator();
        while (e.MoveNext())
        {
            switch (e.Current)
            {
                case "--mum": o = o with { MatchMode = MummerMatchMode.Mum }; break;
                case "--mumcand" or "--mumreference": o = o with { MatchMode = MummerMatchMode.MumReference }; break;
                case "--maxmatch": o = o with { MatchMode = MummerMatchMode.MaxMatch }; break;
                case "-l": o = o with { MinMatch = int.Parse(Next(e), CultureInfo.InvariantCulture) }; break;
                case "-c": o = o with { MinCluster = long.Parse(Next(e), CultureInfo.InvariantCulture) }; break;
                case "-g": o = o with { MaxGap = long.Parse(Next(e), CultureInfo.InvariantCulture) }; break;
                case "-b": o = o with { BreakLength = int.Parse(Next(e), CultureInfo.InvariantCulture) }; break;
                case "-d": o = o with { DiagonalFactor = double.Parse(Next(e), CultureInfo.InvariantCulture) }; break;
                case "-m": o = o with { MaskLength = int.Parse(Next(e), CultureInfo.InvariantCulture) }; break;
                case "-x": o = o with { Matrix = int.Parse(Next(e), CultureInfo.InvariantCulture) }; break;
                case "-p": o = o with { Prefix = Next(e) }; break;
                case "-o": o = o with { Coords = true }; break;
                case "--noextend": o = o with { Extend = false }; break;
                case "--nooptimize": o = o with { Optimize = false }; break;
                case "--nodelta": o = o with { Delta = false }; break;
                default: throw new ArgumentException(e.Current);
            }
        }

        return o;

        static string Next(IEnumerator<string> e) => e.MoveNext() ? e.Current : throw new ArgumentException("missing value");
    }

    private static bool IsReverseAlignmentLine(string line)
    {
        var f = line.Split(' ');
        return f.Length == 7 && long.TryParse(f[2], out long sq) && long.TryParse(f[3], out long eq) && sq > eq;
    }

    public static (List<string> Failures, Dictionary<string, int> Counts) CheckAll(string resource = "promer_oracle.jsonl.gz")
    {
        var failures = new List<string>();
        int rows = 0, errors = 0, withAlignments = 0, reverse = 0, multi = 0, nodelta = 0, coords = 0, matrices = 0, warnings = 0;
        var modes = new HashSet<string>();
        foreach (var j in MummerOracleData.Rows(resource))
        {
            rows++;
            var args = j.StrList("args");
            string label = $"row {rows} [{string.Join(' ', args)}]";
            var o = ParseArgs(args);
            PromerResult? r = null;
            try
            {
                r = PromerProgram.Run(j.Str("ref"), j.Str("qry"), WorkDir + "/ref.fa", WorkDir + "/qry.fa", o);
            }
            catch (Exception e) when (e is ArgumentException or FormatException)
            {
            }

            if (j.GetProperty("rc").GetInt32() != 0)
            {
                errors++;
                if (r is not null)
                    failures.Add(label + ": the program fails, the port returns text");
                continue;
            }

            if (r is null)
            {
                failures.Add(label + ": <error>");
                continue;
            }

            void Compare(string what, string ours, string? expected)
            {
                if (expected is null || ours != expected)
                    failures.Add($"{label}: {what} {MummerOracleData.Diff(ours, expected ?? "<null>")}");
            }

            Compare("aaref", r.AminoReference, j.Str("aaref"));
            Compare("aaqry", r.AminoQuery, j.Str("aaqry"));
            if (j.TryGetProperty("mgaps", out var mg) && mg.ValueKind == JsonValueKind.String)
                Compare("mgaps", r.Mgaps, mg.GetString());
            string? delta = j.StrOrNull("delta"), cluster = j.StrOrNull("cluster"), crd = j.StrOrNull("coords");
            if (delta != r.Delta)
                failures.Add($"{label}: delta {MummerOracleData.Diff(r.Delta ?? "<null>", delta ?? "<null>")}");
            if (cluster != r.Cluster)
                failures.Add($"{label}: cluster {MummerOracleData.Diff(r.Cluster ?? "<null>", cluster ?? "<null>")}");
            if (crd != r.Coords)
                failures.Add($"{label}: coords {MummerOracleData.Diff(r.Coords ?? "<null>", crd ?? "<null>")}");
            if (r.Messages.TrimEnd('\n') != j.Str("err").TrimEnd('\n'))
                failures.Add($"{label}: stderr {MummerOracleData.Diff(r.Messages, j.Str("err"))}");

            modes.Add(o.MatchMode.ToString());
            if (delta != null && delta.Split('\n').Length > 6) withAlignments++;
            if (delta != null)
                reverse += delta.Split('\n').Count(IsReverseAlignmentLine);
            if (j.Str("ref").Split('>').Length > 2 || j.Str("qry").Split('>').Length > 2) multi++;
            if (!o.Delta) nodelta++;
            if (o.Coords) coords++;
            if (o.Matrix != 2) matrices++;
            if (r.Messages.Contains("WARNING", StringComparison.Ordinal)) warnings++;
        }

        var counts = new Dictionary<string, int> { ["rows"] = rows, ["errors"] = errors, ["withAlignments"] = withAlignments, ["reverse"] = reverse, ["multi"] = multi, ["nodelta"] = nodelta, ["coords"] = coords, ["matrices"] = matrices, ["warnings"] = warnings, ["modes"] = modes.Count };
        return (failures, counts);
    }
}
