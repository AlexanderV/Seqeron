using System.IO.Compression;
using System.Text.Json;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// A5-10 (B09 F174): the text of <c>mash bounds</c>, <c>mash dist</c> (<c>-t -v -d -i -n -k -s -S</c>) and <c>mash triangle</c> (<c>-E -v -d -i</c>) from <see cref="MashCommands"/> against the installed Mash 2.3, byte for byte,
/// on 220 random cases per command (TestData/Alignment/mash_commands_oracle.jsonl.gz, scripts/alignment/generate_mash_commands_oracle.py).
/// </summary>
[TestFixture]
public class MashCommands_Oracle_Tests
{
    private static readonly Lazy<List<JsonElement>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.mash_commands_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var rows = new List<JsonElement>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
            rows.Add(JsonDocument.Parse(line).RootElement.Clone());
        return rows;
    });

    private static MashInput[] Files(JsonElement r, string property) => [.. r.GetProperty(property).EnumerateArray().Select(f => new MashInput(f.GetProperty("name").GetString()!,
        [.. f.GetProperty("records").EnumerateArray().Select(x => new MashRecord(x.GetProperty("id").GetString()!, x.GetProperty("comment").GetString()!, x.GetProperty("seq").GetString()!))]))];

    private static double Fraction(JsonElement e) => e.ValueKind == JsonValueKind.Null ? 1.0 : e.GetDouble();

    private static string Produce(JsonElement r)
    {
        string cmd = r.GetProperty("cmd").GetString()!;
        if (cmd == "bounds")
            return MashCommands.Bounds(r.GetProperty("k").GetInt32(), r.GetProperty("p").GetDouble());
        int k = r.GetProperty("k").GetInt32(), s = r.GetProperty("s").GetInt32();
        uint seed = r.GetProperty("seed").GetUInt32();
        bool canonical = !r.GetProperty("noncanonical").GetBoolean(), individual = r.GetProperty("individual").GetBoolean();
        double v = Fraction(r.GetProperty("v")), d = Fraction(r.GetProperty("d"));
        return cmd == "dist"
            ? MashCommands.Dist(Files(r, "refs"), Files(r, "queries"), k, s, seed, canonical, individual, r.GetProperty("table").GetBoolean(), v, d)
            : MashCommands.Triangle(Files(r, "refs"), k, s, seed, canonical, individual, r.GetProperty("edge").GetBoolean(), v, d);
    }

    private static void Check(string command)
    {
        var rows = Oracle.Value.Where(r => r.GetProperty("cmd").GetString() == command).ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(200), command);
        var failures = new List<string>();
        int filtered = 0;
        foreach (var r in rows)
        {
            string expected = r.GetProperty("expected").GetString()!;
            string got = Produce(r);
            if (got != expected)
                failures.Add($"{command} {r.GetRawText()[..Math.Min(120, r.GetRawText().Length)]}\n   expected {expected[..Math.Min(200, expected.Length)]}\n   actual   {got[..Math.Min(200, got.Length)]}");
            if (command != "bounds" && (r.GetProperty("v").ValueKind != JsonValueKind.Null || r.GetProperty("d").ValueKind != JsonValueKind.Null))
                filtered++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(2)));
        if (command != "bounds")
            Assert.That(filtered, Is.GreaterThan(rows.Count / 10), "rows with a -v / -d filter");
    }

    [Test] public void Bounds_EqualsMash_OnAllRows() => Check("bounds");
    [Test] public void Dist_TableAndFilters_EqualsMash_OnAllRows() => Check("dist");
    [Test] public void Triangle_MatrixAndEdgeList_EqualsMash_OnAllRows() => Check("triangle");

    [Test]
    public void Bounds_DefaultParameters_FirstRowsOfTheProgram()
    {
        // mash bounds -k 21 -p 0.99 (2.3)
        string text = MashCommands.Bounds();
        Assert.That(text, Does.StartWith("\nParameters (run with -h for details):\n   k:   21\n   p:   0.99\n\n\tMash distance\nSketch\t0.05\t0.1\t0.15\t0.2\t0.25\t0.3\t0.35\t0.4\n100\t0.0270708\t0.0867606\tinf\tinf"));
        Assert.That(text, Does.Contain("\n1000\t0.00675661\t0.0157591\t0.0322692\t0.0630219\tinf"));
        Assert.Throws<ArgumentOutOfRangeException>(() => MashCommands.Bounds(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => MashCommands.Bounds(21, 1.0));
    }

    [Test]
    public void Dist_ShortSequenceInIndividualMode_IsNotSketched()
    {
        var files = new[] { new MashInput("a.fa", [new MashRecord("long", "", new string('A', 40) + "CGTACGTAGCTAGCTAGCATCGATCGAT"), new MashRecord("short", "", "ACG")]) };
        string text = MashCommands.Dist(files, files, k: 9, sketchSize: 20, individual: true, table: true);
        Assert.That(text, Does.StartWith("#query\tlong\n"));
        Assert.That(text, Does.Not.Contain("short"));
    }
}
