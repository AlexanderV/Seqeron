using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A3-9b (B09 F200 - F202): the delta programs of MUMmer 3.23 on delta files written by nucmer from rearranged, duplicated, inverted, mutated multi-contig inputs (TestData/Alignment/mummer_show_oracle.jsonl.gz,
/// scripts/alignment/generate_mummer_show_oracle.py; the programs ran with <c>time ()</c> fixed to 1700000000 so that <c>srand (time (NULL))</c> and the btab date are reproducible): <c>delta-filter</c>, <c>show-coords</c>,
/// <c>show-aligns</c>, <c>show-snps</c>, <c>show-diff</c>, <c>show-tiling</c> (with its -t / -u / -p files) and <c>dnadiff</c> (all files; the hash-ordered lines of the report and of .unref / .unqry are compared as sets of lines).
/// </summary>
[TestFixture]
public class MummerShow_Oracle_Tests
{
    [Test]
    public void DeltaPrograms_OutputEqualsTheProgram_AllRuns()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.mummer_show_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        var perProgram = new Dictionary<string, int>();
        int cases = 0, errorRuns = 0, tilingWithPath = 0, xmlLinks = 0, lasRuns = 0, snpSelect = 0, filterLis = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
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

        Assert.That(failures, Is.Empty, $"{failures.Count} runs differ:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(cases, Is.GreaterThanOrEqualTo(150));
        foreach (string program in new[] { "delta-filter", "show-coords", "show-aligns", "show-snps", "show-diff", "show-tiling", "dnadiff" })
            Assert.That(perProgram.GetValueOrDefault(program), Is.GreaterThanOrEqualTo(200), program);
        Assert.That(tilingWithPath, Is.GreaterThan(100));
        Assert.That(xmlLinks, Is.GreaterThan(3));
        Assert.That(lasRuns, Is.GreaterThan(3));
        Assert.That(snpSelect, Is.GreaterThan(5));
        Assert.That(filterLis, Is.GreaterThan(100));
    }
}
