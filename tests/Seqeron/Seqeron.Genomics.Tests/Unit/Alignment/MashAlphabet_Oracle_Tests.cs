using System.Globalization;
using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-10b2 (B09 F214): the alphabet options of Mash 2.3 - <c>-a</c> (amino acids, implies <c>-n -k 9</c>), <c>-z TEXT</c> (an own alphabet, implies <c>-n</c>), <c>-Z</c> (the case is kept) - for <c>mash sketch</c> (the whole .msh file),
/// <c>dist</c> / <c>triangle</c> (the text) and <c>screen</c> (a protein query sketch translates the mixture in six frames) against the installed mash 2.3
/// (TestData/Alignment/mash_alphabet_oracle.jsonl.gz, scripts/alignment/generate_mash_alphabet_oracle.py).
/// </summary>
[TestFixture]
public class MashAlphabet_Oracle_Tests
{
    private sealed record Row(int Id, string Kind, List<string> Args, List<(string Name, string Fasta)> Files, byte[]? QueryMsh, string Out, byte[]? Msh);

    private static readonly Lazy<List<Row>> Oracle = new(Load);

    private static List<Row> Load()
    {
        var rows = new List<Row>();
        foreach (string part in new[] { "sketch", "dist", "screen" })
            Load(part, rows);
        return rows;
    }

    private static void Load(string part, List<Row> rows)
    {
        string name = $"Seqeron.Genomics.Tests.TestData.Alignment.mash_alphabet_{part}_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            string? q = j.GetProperty("qmsh").ValueKind == JsonValueKind.Null ? null : j.GetProperty("qmsh").GetString();
            string? m = j.GetProperty("msh").ValueKind == JsonValueKind.Null ? null : j.GetProperty("msh").GetString();
            rows.Add(new Row(j.GetProperty("id").GetInt32(), j.GetProperty("kind").GetString()!, [.. j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!)],
                [.. j.GetProperty("files").EnumerateArray().Select(f => (f.GetProperty("name").GetString()!, f.GetProperty("fasta").GetString()!))],
                q is null ? null : Convert.FromBase64String(q), j.GetProperty("out").GetString()!, m is null ? null : Convert.FromBase64String(m)));
        }
    }

    // the command line of a row: the value options (the screen's -i is a number, the others' a flag), the flags and the files
    private sealed class Cli
    {
        public readonly Dictionary<string, string> Values = [];
        public readonly HashSet<string> Flags = [];
        public readonly List<string> Positional = [];

        public Cli(IReadOnlyList<string> args, bool iTakesValue)
        {
            var valued = new HashSet<string> { "-k", "-s", "-S", "-m", "-c", "-g", "-b", "-I", "-z", "-v", "-d", "-o" };
            if (iTakesValue)
                valued.Add("-i");
            var queue = new Queue<string>(args.Skip(1));
            while (queue.Count > 0)
            {
                string a = queue.Dequeue();
                if (valued.Contains(a))
                    Values[a] = queue.Dequeue();
                else if (a.StartsWith('-') && a.Length == 2)
                    Flags.Add(a);
                else
                    Positional.Add(a);
            }
        }

        public int Int(string flag, int fallback) => Values.TryGetValue(flag, out var v) ? int.Parse(v, CultureInfo.InvariantCulture) : fallback;
    }

    private static MashInput Input(Row r, string name) => MashCommands.ParseFasta(name, r.Files.First(f => f.Name == name).Fasta);

    private static (int K, MashAlphabet Alphabet) KAndAlphabet(Cli c)
    {
        var alphabet = MashAlphabet.FromOptions(c.Flags.Contains("-a"), c.Values.GetValueOrDefault("-z"), c.Flags.Contains("-Z"));
        return (MashAlphabet.ResolveKmerSize(c.Values.TryGetValue("-k", out var k) ? int.Parse(k, CultureInfo.InvariantCulture) : null, c.Flags.Contains("-a")), alphabet);
    }

    private static MashCommands.MashReadsOptions? Reads(Cli c)
    {
        long bloom = c.Values.TryGetValue("-b", out var b) ? MashCommands.ParseSize(b) : 0;
        if (!c.Flags.Contains("-r") && !c.Values.ContainsKey("-m") && !c.Values.ContainsKey("-c") && !c.Values.ContainsKey("-g") && bloom == 0)
            return null;
        return new MashCommands.MashReadsOptions(c.Int("-m", 1), c.Values.TryGetValue("-c", out var cov) ? double.Parse(cov, CultureInfo.InvariantCulture) : 0, c.Values.TryGetValue("-g", out var g) ? long.Parse(g, CultureInfo.InvariantCulture) : 0, bloom);
    }

    private static string Run(Row r)
    {
        var c = new Cli(r.Args, r.Kind == "screen");
        if (r.Kind == "screen")
        {
            var queryFile = MashSketchFile.Read(r.QueryMsh!);
            var mixtures = c.Positional.Skip(1).Select(n => Input(r, n)).ToList();
            double identity = c.Values.TryGetValue("-i", out var i) ? double.Parse(i, CultureInfo.InvariantCulture) : 0.0;
            double pvalue = c.Values.TryGetValue("-v", out var v) ? double.Parse(v, CultureInfo.InvariantCulture) : 1.0;
            return MashCommands.Screen([new MashInput("q.msh", [], queryFile)], mixtures, winnerTakesAll: c.Flags.Contains("-w"), minIdentity: identity, maxPValue: pvalue);
        }

        var (k, alphabet) = KAndAlphabet(c);
        int s = c.Int("-s", 1000);
        uint seed = MashCommands.SeedFromOption(c.Values.TryGetValue("-S", out var sd) ? long.Parse(sd, CultureInfo.InvariantCulture) : 42);
        bool canonical = !c.Flags.Contains("-n");
        bool individual = c.Flags.Contains("-i");
        var reads = Reads(c);
        switch (r.Kind)
        {
            case "dist":
                {
                    var inputs = c.Positional.Select(n => Input(r, n)).ToList();
                    return MashCommands.Dist([inputs[0]], inputs.Skip(1).ToList(), k, s, seed, canonical, individual, c.Flags.Contains("-t"),
                        c.Values.TryGetValue("-v", out var pv) ? double.Parse(pv, CultureInfo.InvariantCulture) : 1.0, c.Values.TryGetValue("-d", out var dd) ? double.Parse(dd, CultureInfo.InvariantCulture) : 1.0,
                        default, reads, c.Flags.Contains("-C"), alphabet);
                }
            case "triangle":
                return MashCommands.Triangle(c.Positional.Select(n => Input(r, n)).ToList(), k, s, seed, canonical, individual, c.Flags.Contains("-E"),
                    c.Values.TryGetValue("-v", out var tv) ? double.Parse(tv, CultureInfo.InvariantCulture) : 1.0, c.Values.TryGetValue("-d", out var td) ? double.Parse(td, CultureInfo.InvariantCulture) : 1.0,
                    default, reads, c.Flags.Contains("-C"), alphabet);
            default:
                throw new InvalidOperationException(r.Kind);
        }
    }

    private static byte[] Sketch(Row r, out MashSketchFileData data)
    {
        var c = new Cli(r.Args, false);
        var (k, alphabet) = KAndAlphabet(c);
        var inputs = c.Positional.Select(n => Input(r, n)).ToList();
        data = MashCommands.Sketch(inputs, k, c.Int("-s", 1000), MashCommands.SeedFromOption(c.Values.TryGetValue("-S", out var sd) ? long.Parse(sd, CultureInfo.InvariantCulture) : 42), !c.Flags.Contains("-n"), c.Flags.Contains("-i"), default,
            Reads(c), c.Values.GetValueOrDefault("-I"), null, alphabet);
        return MashSketchFile.Write(data);
    }

    [Test]
    public void SketchFile_WithAlphabetOptions_EqualsMash_ByteForByte()
    {
        var rows = Oracle.Value.Where(r => r.Kind == "sketch").ToList();
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(200));
        var failures = new List<string>();
        int protein = 0, custom = 0, preserve = 0, wide = 0, readsRows = 0;
        foreach (var r in rows)
        {
            byte[] got = Sketch(r, out var data);
            if (!got.AsSpan().SequenceEqual(r.Msh!))
                failures.Add($"row {r.Id} {string.Join(' ', r.Args)}: expected {r.Msh!.Length} bytes, got {got.Length}");
            else if (MashCommands.Info(data) != r.Out)
                failures.Add($"row {r.Id}: mash info differs");
            if (r.Args.Contains("-a")) protein++;
            if (r.Args.Contains("-z")) custom++;
            if (r.Args.Contains("-Z")) preserve++;
            if (data.Use64) wide++;
            if (r.Args.Contains("-r") || r.Args.Contains("-b")) readsRows++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(4)));
        Assert.That(protein, Is.GreaterThanOrEqualTo(60));
        Assert.That(custom, Is.GreaterThanOrEqualTo(60));
        Assert.That(preserve, Is.GreaterThanOrEqualTo(40));
        Assert.That(wide, Is.GreaterThanOrEqualTo(40), "64-bit sketches");
        Assert.That(readsRows, Is.GreaterThanOrEqualTo(30), "read-set sketches");
    }

    [Test]
    public void DistAndTriangle_WithAlphabetOptions_EqualMash()
    {
        var rows = Oracle.Value.Where(r => r.Kind is "dist" or "triangle").ToList();
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(200));
        var failures = new List<string>();
        int withOutput = 0;
        foreach (var r in rows)
        {
            string got = Run(r);
            if (got != r.Out)
                failures.Add($"row {r.Id} {string.Join(' ', r.Args)}\n--- expected\n{r.Out}--- actual\n{got}");
            if (r.Out.Length > 0)
                withOutput++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(2)));
        Assert.That(withOutput, Is.GreaterThanOrEqualTo(rows.Count * 3 / 4));
    }

    [Test]
    public void Screen_WithAlphabetQueries_EqualsMash()
    {
        var rows = Oracle.Value.Where(r => r.Kind == "screen").ToList();
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(200));
        var failures = new List<string>();
        int protein = 0, withOutput = 0;
        foreach (var r in rows)
        {
            string got = Run(r);
            if (got != r.Out)
                failures.Add($"row {r.Id} {string.Join(' ', r.Args)}\n--- expected\n{r.Out}--- actual\n{got}");
            if (MashSketchFile.Read(r.QueryMsh!).Alphabet == MashAlphabet.ProteinCharacters)
                protein++;
            if (r.Out.Length > 0)
                withOutput++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(2)));
        Assert.That(protein, Is.GreaterThanOrEqualTo(150), "translated (amino-acid) screens");
        Assert.That(withOutput, Is.GreaterThanOrEqualTo(rows.Count / 2));
    }
}
