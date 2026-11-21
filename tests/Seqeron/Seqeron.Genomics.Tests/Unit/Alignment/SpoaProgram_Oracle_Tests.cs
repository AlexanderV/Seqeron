using System.Text;
using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-26c, B09 F232): the <c>spoa</c> 4.0.8 command line (<see cref="SpoaProgram"/>) equals the installed Debian program run on the same arguments and files - exit status (134 for an abort),
/// standard output, standard error and the written DOT files (scripts/alignment/generate_spoa_oracle.py, spoa_oracle.tsv).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SpoaProgram_Oracle_Tests
{
    private sealed record Row(string Id, List<string> Args, Dictionary<string, byte[]> Files, int Rc, string StdOut, string StdErr, Dictionary<string, string> OutFiles);

    private static string Unesc(string s) => System.Text.RegularExpressions.Regex.Replace(s, @"\\(.)", m => m.Groups[1].Value switch { "t" => "\t", "n" => "\n", "r" => "\r", var c => c });

    private static readonly Lazy<List<Row>> Rows = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.spoa_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            var args = JsonNode.Parse(Unesc(f[1]))!.AsArray().Select(x => x!.GetValue<string>()).ToList();
            var files = JsonNode.Parse(Unesc(f[2]))!.AsObject().ToDictionary(kv => kv.Key, kv => Convert.FromBase64String(kv.Value!.GetValue<string>()));
            var outFiles = JsonNode.Parse(Unesc(f[6]))!.AsObject().ToDictionary(kv => kv.Key, kv => kv.Value!.GetValue<string>());
            list.Add(new Row(f[0], args, files, int.Parse(f[3]), Unesc(f[4]), Unesc(f[5]), outFiles));
        }

        return list;
    });

    [Test]
    public void Oracle_HasAtLeast200Runs()
    {
        Assert.That(Rows.Value.Count, Is.GreaterThanOrEqualTo(200));
    }

    [Test]
    public void Run_EveryCommandLine_EqualsTheInstalledProgram() => CheckAll(256);

    [Test]
    public void Run_With128BitVectors_GivesTheSameResult() => CheckAll(128);

    private static void CheckAll(int registerBits)
    {
        var failures = new List<string>();
        foreach (var r in Rows.Value)
        {
            var result = SpoaProgram.Run(r.Args, p => r.Files.TryGetValue(p, out var b) ? b : null, registerBits);
            var written = result.Files.ToDictionary(kv => kv.Key, kv => kv.Value);
            if (result.ExitCode != r.Rc || result.StdOut != r.StdOut || result.StdErr != r.StdErr || !written.Count.Equals(r.OutFiles.Count) || r.OutFiles.Any(kv => !written.TryGetValue(kv.Key, out var v) || v != kv.Value))
                failures.Add($"{r.Id} spoa {string.Join(' ', r.Args)}: rc {result.ExitCode} / {r.Rc}\n--- expected stdout\n{Trim(r.StdOut)}\n--- actual stdout\n{Trim(result.StdOut)}\n--- expected stderr {Trim(r.StdErr)}\n--- actual stderr {Trim(result.StdErr)}");
        }

        Assert.That(failures, Is.Empty, string.Join("\n=====\n", failures.Take(3)) + $"\n({failures.Count} of {Rows.Value.Count}: {string.Join(", ", failures.Take(60).Select(f => f.Split(' ')[0]))})");
    }

    private static string Trim(string s) => s.Length > 1500 ? s[..1500] + "..." : s;
}
