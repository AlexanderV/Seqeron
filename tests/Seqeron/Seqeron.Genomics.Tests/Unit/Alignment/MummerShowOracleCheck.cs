using System.Text.Json;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>The oracle comparison of the MUMmer delta programs on a recorded set of delta files (no test framework: also run from a scratch console).</summary>
public static class MummerShowOracleCheck
{
    public static (List<string> Failures, int Cases, Dictionary<string, int> PerProgram, Dictionary<string, int> Features) CheckAll(string resource)
    {
        var failures = new List<string>();
        var perProgram = new Dictionary<string, int>();
        int cases = 0, errorRuns = 0, tilingWithPath = 0, xmlLinks = 0, lasRuns = 0, snpSelect = 0, filterLis = 0;
        foreach (var j in MummerOracleData.Rows(resource))
        {
            cases++;
            string refFasta = j.GetProperty("ref").GetString()!, qryFasta = j.GetProperty("qry").GetString()!, delta = j.GetProperty("delta").GetString()!;
            int runNo = 0;
            foreach (var r in j.GetProperty("runs").EnumerateArray())
            {
                runNo++;
                string prog = r.GetProperty("prog").GetString()!;
                var args = r.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
                string? stdin = r.GetProperty("stdin").ValueKind == JsonValueKind.Null ? null : r.GetProperty("stdin").GetString();
                string expected = r.GetProperty("out").GetString()!;
                int rc = r.GetProperty("rc").GetInt32();
                string? got = MummerShowOracleRunner.Run(prog, args, stdin, delta, refFasta, qryFasta, out var files);
                string label = $"case {cases} run {runNo} {prog} {string.Join(' ', args)}";
                perProgram[prog] = perProgram.GetValueOrDefault(prog) + 1;
                if (rc != 0)
                {
                    errorRuns++;
                    if (got is not null)
                        failures.Add(label + ": the program fails, the port returns text");
                    continue;
                }

                if (got != expected)
                {
                    failures.Add($"{label}: {(got is null ? "<error>" : MummerShowOracleRunner.Diff(got, expected))}");
                    continue;
                }

                if (r.TryGetProperty("files", out var ef) && ef.ValueKind == JsonValueKind.Object)
                {
                    foreach (var f in ef.EnumerateObject())
                    {
                        string? ex = f.Value.ValueKind == JsonValueKind.Null ? null : f.Value.GetString();
                        bool unordered = prog == "dnadiff" && f.Name is "report" or "unref" or "unqry";
                        if (!MummerShowOracleRunner.FileEqual(f.Name, files.GetValueOrDefault(f.Name), ex, unordered))
                            failures.Add($"{label}: file {f.Name} differs: {MummerShowOracleRunner.Diff(files.GetValueOrDefault(f.Name) ?? "<null>", ex ?? "<null>")}");
                    }
                }

                if (prog == "show-tiling" && expected.Contains('\t')) tilingWithPath++;
                if (prog == "show-tiling" && expected.Contains("<LINK", StringComparison.Ordinal)) xmlLinks++;
                if (prog == "show-coords" && (args.Contains("-g") || args.Contains("-G"))) lasRuns++;
                if (prog == "show-snps" && args.Contains("-S")) snpSelect++;
                if (prog == "delta-filter" && args.Any(a => a is "-q" or "-r" or "-g" or "-1" or "-m")) filterLis++;
            }
        }


        var features = new Dictionary<string, int> { ["errorRuns"] = errorRuns, ["tilingWithPath"] = tilingWithPath, ["xmlLinks"] = xmlLinks, ["lasRuns"] = lasRuns, ["snpSelect"] = snpSelect, ["filterLis"] = filterLis };
        return (failures, cases, perProgram, features);
    }
}
