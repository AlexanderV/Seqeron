using System.IO.Compression;
using System.Reflection;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// B09 WP61 (F280): the EMBOSS 6.6.0 sequence programs run from their command line (<see cref="EmbossProgramRunner"/>) equal the installed programs byte for byte: the standard output, the error stream, the exit status and
/// every file the program writes (with the clock of the recorded run: <see cref="EmbossRunEnvironment"/>), on the rows of TestData/Alignment/emboss_prog_&lt;program&gt;_oracle.jsonl.gz (scripts/alignment/emboss_prog_oracle_lib.py and the generate_emboss_prog_*_oracle.py scripts of the programs:
/// every row is one run of <c>&lt;program&gt; -auto ARGS</c> in a work directory that holds the input files). A row on which the reference program dies by a signal or never ends must be refused by the port
/// (<see cref="InvalidOperationException"/>: the result is not defined).
/// </summary>
[TestFixture]
public class EmbossPrograms_Oracle_Tests
{
    private const string Prefix = "emboss_prog_";
    private const string Suffix = "_oracle.jsonl.gz";

    public static IEnumerable<string> Programs()
    {
        foreach (string name in Assembly.GetExecutingAssembly().GetManifestResourceNames().OrderBy(n => n, StringComparer.Ordinal))
        {
            int i = name.IndexOf(Prefix, StringComparison.Ordinal);
            if (i >= 0 && name.EndsWith(Suffix, StringComparison.Ordinal))
                yield return name[(i + Prefix.Length)..^Suffix.Length];
        }
    }

    [TestCaseSource(nameof(Programs))]
    public void Program_EqualsTheInstalledProgram_AllRows(string program)
    {
        string resource = $"Seqeron.Genomics.Tests.TestData.Alignment.{Prefix}{program}{Suffix}";
        using var stream = Assembly.GetExecutingAssembly().GetManifestResourceStream(resource) ?? throw new InvalidOperationException($"Embedded oracle '{resource}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, nonEmpty = 0, died = 0, refused = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            var files = new Dictionary<string, string>();
            foreach (var p in j.GetProperty("files").EnumerateObject())
                files[p.Name] = p.Value.GetString()!;
            int rc = j.GetProperty("rc").GetInt32();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc < 0;
            string id = $"row {j.GetProperty("id")} [{string.Join(' ', args)}]";
            // the clock of the recorded run: the date the writers print ('epoch', UTC) and, for the random programs, the microseconds of gettimeofday ('usec') or the variable timetoday
            var clock = new EmbossRunEnvironment(
                j.TryGetProperty("epoch", out var epoch) ? DateTimeOffset.FromUnixTimeSeconds(epoch.GetInt64()).UtcDateTime : null,
                j.TryGetProperty("usec", out var usec) && usec.ValueKind == JsonValueKind.Number ? usec.GetInt32() : null,
                j.TryGetProperty("timetoday", out var timetoday) && timetoday.ValueKind == JsonValueKind.String ? timetoday.GetString() : null)
            {
                // the run had its data directory WORK/edata/ (EMBOSS_DATA; the PROSITE programs: scripts/alignment/emboss_prog_oracle_lib.py, run(edata=True))
                DataDirectory = j.TryGetProperty("edata", out var edata) && edata.ValueKind == JsonValueKind.True ? "edata/" : null,
            };
            try
            {
                var r = EmbossProgramRunner.Run(program, args, files, clock);
                if (refusal)
                {
                    failures.Add($"{id}: the reference dies by a signal or never ends (status {rc}), the port returned");
                    continue;
                }

                var expectedFiles = new SortedDictionary<string, string>(StringComparer.Ordinal);
                foreach (var p in j.GetProperty("outfiles").EnumerateObject())
                    expectedFiles[p.Name] = p.Value.GetString()!;
                string why = "";
                if (r.ExitCode != rc)
                    why = $"exit status {r.ExitCode} expected {rc}";
                else if (r.Stdout != j.GetProperty("out").GetString())
                    why = "standard output differs";
                else if (r.Stderr != j.GetProperty("err").GetString())
                    why = $"error stream '{r.Stderr.Replace("\n", "|")}' expected '{j.GetProperty("err").GetString()!.Replace("\n", "|")}'";
                else if (!expectedFiles.SequenceEqual(new SortedDictionary<string, string>(r.Files.ToDictionary(k => k.Key, k => k.Value), StringComparer.Ordinal)))
                    why = $"files {string.Join(",", r.Files.Keys)} expected {string.Join(",", expectedFiles.Keys)} (or their text differs)";
                if (why.Length > 0)
                    failures.Add($"{id}: {why}");
                else
                {
                    if (r.Stdout.Length > 0 || r.Files.Count > 0)
                        nonEmpty++;
                    if (rc != 0)
                        died++;
                }
            }
            catch (InvalidOperationException) when (refusal)
            {
                refused++;
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(12)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(200), "at least 200 runs of every program");
        Assert.That(nonEmpty, Is.GreaterThan(rows / 4), "runs that wrote output");
        _ = died;
        _ = refused;
    }
}
