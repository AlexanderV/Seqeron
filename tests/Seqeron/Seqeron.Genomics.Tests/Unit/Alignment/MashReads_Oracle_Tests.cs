using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-10b (B09 F186): the read-set mode of Mash 2.3 <c>sketch</c> / <c>dist</c> / <c>triangle</c> (<c>-r</c> with <c>-m</c> minimum copies = the noise filter, <c>-c</c> target coverage, <c>-g</c> genome size; <see cref="MashCommands.MashReadsOptions"/>)
/// and <c>-C</c> (comments with the names), as TEXT against the installed mash 2.3 on 260 random runs (TestData/Alignment/mash_reads_oracle.jsonl.gz, scripts/alignment/generate_mash_reads_oracle.py: read sets of 1 - 12x coverage with
/// errors, k 4 - 32, sketch sizes 5 - 2000, seeds, <c>-n</c>, <c>-v</c>, <c>-d</c>, <c>-t</c>, <c>-E</c>). The Bloom filter <c>-b</c> is not part of the port.
/// </summary>
[TestFixture]
public class MashReads_Oracle_Tests
{
    private static List<MashInput> Inputs(JsonElement array) => [.. array.EnumerateArray().Select(f => new MashInput(f.GetProperty("name").GetString()!,
        [.. f.GetProperty("records").EnumerateArray().Select(r => new MashRecord(r.GetProperty("id").GetString()!, r.GetProperty("comment").GetString()!, r.GetProperty("seq").GetString()!))]))];

    [Test]
    public void ReadsMode_AndComments_EqualMash()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.mash_reads_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        int rows = 0, readRows = 0, filtered = 0, comments = 0, covered = 0;
        var failures = new List<string>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            string cmd = j.GetProperty("cmd").GetString()!;
            int k = j.GetProperty("k").GetInt32(), s = j.GetProperty("s").GetInt32();
            uint seed = (uint)j.GetProperty("seed").GetInt64();
            bool canonical = !j.GetProperty("noncanonical").GetBoolean();
            MashCommands.MashReadsOptions? reads = j.GetProperty("reads").GetBoolean()
                ? new MashCommands.MashReadsOptions(j.GetProperty("minCopies").GetInt32(), j.GetProperty("targetCoverage").GetDouble(), j.GetProperty("genomeSize").GetInt64()) : null;
            bool showComments = j.GetProperty("comments").GetBoolean();
            double v = j.GetProperty("v").GetDouble(), d = j.GetProperty("d").GetDouble();
            string got = cmd == "dist"
                ? MashCommands.Dist(Inputs(j.GetProperty("refs")), Inputs(j.GetProperty("queries")), k, s, seed, canonical, false, j.GetProperty("table").GetBoolean(), v, d, default, reads, showComments)
                : MashCommands.Triangle(Inputs(j.GetProperty("refs")), k, s, seed, canonical, false, j.GetProperty("edge").GetBoolean(), v, d, default, reads, showComments);
            string expected = j.GetProperty("expected").GetString()!;
            if (got != expected)
                failures.Add($"row {rows} {cmd} k={k} s={s}\n--- expected\n{expected}--- actual\n{got}");
            rows++;
            if (reads is not null)
                readRows++;
            if (reads is { MinCopies: > 1 })
                filtered++;
            if (reads is { TargetCoverage: > 0 })
                covered++;
            if (showComments)
                comments++;
        }

        Assert.That(rows, Is.GreaterThanOrEqualTo(200));
        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(2)));
        Assert.That((readRows, filtered, covered, comments), Is.GreaterThanOrEqualTo((150, 60, 40, 40)).Using<(int, int, int, int)>((a, b) => a.Item1 >= b.Item1 && a.Item2 >= b.Item2 && a.Item3 >= b.Item3 && a.Item4 >= b.Item4 ? 1 : -1), "rows per option");
    }

    [Test]
    public void ReadsMode_IsRefusedWithIndividualSketches()
    {
        var f = new MashInput("a.fa", [new MashRecord("a", "", "ACGTACGTACGTAGCTAGCTAGGATCCA")]);
        Assert.Throws<ArgumentException>(() => MashCommands.Dist([f], [f], 8, 20, 42, true, true, false, 1.0, 1.0, default, new MashCommands.MashReadsOptions()));
    }
}
