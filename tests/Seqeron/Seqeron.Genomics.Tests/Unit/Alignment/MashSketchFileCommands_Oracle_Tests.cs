using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-10b (B09 F193): <c>mash paste</c> and <c>mash dist</c> / <c>triangle</c> / <c>screen</c> with sketch files (.msh) among their inputs, against the installed mash 2.3 on 160 random runs
/// (TestData/Alignment/mash_msh_commands_oracle.jsonl.gz, scripts/alignment/generate_mash_msh_commands_oracle.py): the parameters of the first sketch file are the command's, a later sketch file with another
/// seed / k / strand mode or a smaller sketch size is skipped and a larger one reduced, sequence files are sketched with the inherited parameters; the pasted file is compared byte for byte.
/// </summary>
[TestFixture]
public class MashSketchFileCommands_Oracle_Tests
{
    private sealed record Row(string Cmd, List<string> Args, List<(string Name, string Kind, string Content)> Inputs, string? Expected, byte[]? Paste);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.mash_msh_commands_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows.Add(new Row(j.GetProperty("cmd").GetString()!, [.. j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!)],
                [.. j.GetProperty("inputs").EnumerateArray().Select(i => (i.GetProperty("name").GetString()!, i.GetProperty("kind").GetString()!, i.GetProperty("kind").GetString() == "msh" ? i.GetProperty("msh").GetString()! : i.GetProperty("fasta").GetString()!))],
                j.TryGetProperty("expected", out var e) ? e.GetString() : null, j.TryGetProperty("paste", out var p) ? Convert.FromBase64String(p.GetString()!) : null));
        }

        return rows;
    });

    private static List<MashInput> Inputs(Row r) => [.. r.Inputs.Select(i => i.Kind == "msh"
        ? new MashInput(i.Name, [], MashSketchFile.Read(Convert.FromBase64String(i.Content)))
        : MashSketchFile_Oracle_Tests.ToInput(i.Name, i.Content))];

    private static (int K, int S, uint Seed, bool Canonical) Options(Row r)
    {
        int k = 21, s = 1000;
        uint seed = 42;
        string Value(string flag) => r.Args[r.Args.IndexOf(flag) + 1];
        if (r.Args.Contains("-k")) k = int.Parse(Value("-k"));
        if (r.Args.Contains("-s")) s = int.Parse(Value("-s"));
        if (r.Args.Contains("-S")) seed = MashCommands.SeedFromOption(long.Parse(Value("-S")));
        return (k, s, seed, !r.Args.Contains("-n"));
    }

    [Test]
    public void SketchFileInputs_EqualMash()
    {
        var failures = new List<string>();
        int dist = 0, triangle = 0, screen = 0, paste = 0, mixed = 0;
        foreach (var r in Oracle.Value)
        {
            var (k, s, seed, canonical) = Options(r);
            var inputs = Inputs(r);
            string got;
            if (r.Cmd == "paste")
            {
                paste++;
                var bytes = MashSketchFile.Write(MashCommands.Paste([.. inputs.Select(i => i.SketchFile!)]));
                if (!bytes.AsSpan().SequenceEqual(r.Paste))
                    failures.Add($"paste of {string.Join(' ', r.Inputs.Select(i => i.Name))}: {bytes.Length} vs {r.Paste!.Length} bytes");
                continue;
            }

            if (r.Cmd == "dist")
            {
                dist++;
                got = MashCommands.Dist(inputs.Take(1).ToList(), inputs.Skip(1).ToList(), k, s, seed, canonical, false, r.Args.Contains("-t"));
            }
            else if (r.Cmd == "triangle")
            {
                triangle++;
                got = MashCommands.Triangle(inputs, k, s, seed, canonical, false, r.Args.Contains("-E"));
            }
            else
            {
                screen++;
                got = MashCommands.Screen(inputs.Take(1).ToList(), inputs.Skip(1).ToList(), k, s, seed, canonical);
            }

            if (r.Inputs.Any(i => i.Kind == "seq") && r.Inputs.Any(i => i.Kind == "msh"))
                mixed++;
            if (got != r.Expected)
                failures.Add($"{r.Cmd} {string.Join(' ', r.Args)}\n--- expected\n{r.Expected}--- actual\n{got}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count}:\n" + string.Join("\n", failures.Take(2)));
        Assert.That(dist, Is.GreaterThanOrEqualTo(40));
        Assert.That(triangle, Is.GreaterThanOrEqualTo(30));
        Assert.That(screen, Is.GreaterThanOrEqualTo(30));
        Assert.That(paste, Is.GreaterThanOrEqualTo(30));
        Assert.That(mixed, Is.GreaterThanOrEqualTo(30));
    }
}
