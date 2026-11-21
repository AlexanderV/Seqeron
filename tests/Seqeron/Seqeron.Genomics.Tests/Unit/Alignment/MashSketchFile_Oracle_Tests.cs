using System.Globalization;
using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-10b (B09 F193): the Mash 2.3 sketch FILE. <c>mash sketch -o x.msh</c> writes a Cap'n Proto message of <c>MinHash.capnp</c> built by <c>capnp::MallocMessageBuilder</c>; <see cref="MashSketchFile"/> reproduces its
/// allocation (segments of 1 024, 1 024, 2 048, ... words, far pointers with landing pads, the allocation order of <c>Sketch::writeToCapnp</c>) byte for byte on 130 random runs (genomes and read sets, k 4 - 32, sketch sizes
/// 1 - 12 000, seeds, <c>-n</c>, <c>-i</c>, <c>-r -m -c -g</c>, <c>-l</c>: TestData/Alignment/mash_msh_oracle.jsonl.gz, scripts/alignment/generate_mash_msh_oracle.py) and the text of <c>mash info</c> (default, <c>-t</c>, <c>-H</c>,
/// <c>-d</c>, <c>-c</c>) on the file.
/// </summary>
[TestFixture]
public class MashSketchFile_Oracle_Tests
{
    private sealed record Row(int Id, List<string> Args, List<(string Name, string Fasta)> Files, byte[] Msh, Dictionary<string, string> Info);

    private static readonly Lazy<List<Row>> Oracle = new(() => Load("mash_msh_oracle"));

    private static readonly Lazy<List<Row>> ReadsOracle = new(() => Load("mash_reads_msh_oracle"));

    private static List<Row> Load(string resource)
    {
        string name = $"Seqeron.Genomics.Tests.TestData.Alignment.{resource}.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows.Add(new Row(j.GetProperty("id").GetInt32(), [.. j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!)],
                [.. j.GetProperty("files").EnumerateArray().Select(f => (f.GetProperty("name").GetString()!, f.GetProperty("fasta").GetString()!))],
                Convert.FromBase64String(j.GetProperty("msh").GetString()!), j.GetProperty("info").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!)));
        }

        return rows;
    }

    // FASTA text -> records as Mash reads them (kseq, including the stale comment of a record without one)
    internal static MashInput ToInput(string name, string fasta) => MashCommands.ParseFasta(name, fasta);

    private static byte[] Make(Row r, out MashSketchFileData data)
    {
        int k = 21, s = 1000;
        uint seed = 42;
        bool canonical = true, individual = false;
        MashCommands.MashReadsOptions? reads = null;
        int minCopies = 1;
        double coverage = 0;
        long genome = 0;
        bool readMode = false;
        string Value(string flag) => r.Args[r.Args.IndexOf(flag) + 1];
        if (r.Args.Contains("-k")) k = int.Parse(Value("-k"), CultureInfo.InvariantCulture);
        if (r.Args.Contains("-s")) s = int.Parse(Value("-s"), CultureInfo.InvariantCulture);
        if (r.Args.Contains("-S")) seed = MashCommands.SeedFromOption(long.Parse(Value("-S"), CultureInfo.InvariantCulture));
        if (r.Args.Contains("-m")) minCopies = int.Parse(Value("-m"), CultureInfo.InvariantCulture);
        if (r.Args.Contains("-c")) coverage = double.Parse(Value("-c"), CultureInfo.InvariantCulture);
        if (r.Args.Contains("-g")) genome = long.Parse(Value("-g"), CultureInfo.InvariantCulture);
        canonical = !r.Args.Contains("-n");
        individual = r.Args.Contains("-i");
        readMode = r.Args.Contains("-r");
        long bloom = r.Args.Contains("-b") ? MashCommands.ParseSize(Value("-b")) : 0;

        if (readMode || bloom > 0)
            reads = new MashCommands.MashReadsOptions(minCopies, coverage, genome, bloom);
        var inputs = r.Files.Select(f => ToInput(f.Name, f.Fasta)).ToList();
        string? id = r.Args.Contains("-I") ? Value("-I") : null, comment = r.Args.Contains("-C") ? Value("-C") : null;
        data = MashCommands.Sketch(inputs, k, s, seed, canonical, individual, default, reads, id, comment);
        return MashSketchFile.Write(data);
    }

    /// <summary>A5-10b3 (F194): 140 read sets (10 - 1 200 reads, k 4 - 32, sketch sizes 1 - 3 000, -m 1 - 3, -c, -g, -n, 1 - 2 files): the whole .msh file, the order of the multiplicities (robin_hood iteration order) included.</summary>
    [Test]
    public void ReadSetSketchFile_EqualsMash_ByteForByte_CountsInRobinHoodOrder()
    {
        var rows = ReadsOracle.Value;
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(100));
        var failures = new List<string>();
        foreach (var r in rows)
        {
            var got = Make(r, out var data);
            if (!got.AsSpan().SequenceEqual(r.Msh))
                failures.Add($"row {r.Id} {string.Join(' ', r.Args)}: expected {r.Msh.Length} bytes, got {got.Length}; first difference at {FirstDifference(got, r.Msh)}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(5)));
    }

    private static readonly Lazy<List<Row>> BloomOracle = new(() => Load("mash_bloom_msh_oracle"));

    /// <summary>A5-10b2 (F197): 130 read sets filtered by a Bloom filter (<c>mash sketch -b SIZE</c>, 8 bytes - 2 MB, with and without <c>-r</c>): the whole .msh file byte for byte.</summary>
    [Test]
    public void BloomFilterSketchFile_EqualsMash_ByteForByte()
    {
        var rows = BloomOracle.Value;
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(100));
        var failures = new List<string>();
        foreach (var r in rows)
        {
            var got = Make(r, out _);
            if (!got.AsSpan().SequenceEqual(r.Msh))
                failures.Add($"row {r.Id} {string.Join(' ', r.Args)}: expected {r.Msh.Length} bytes, got {got.Length}; first difference at {FirstDifference(got, r.Msh)}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(5)));
    }

    [Test]
    public void ParseSize_AsMashCommandCpp()
    {
        Assert.That(MashCommands.ParseSize("500"), Is.EqualTo(500));
        Assert.That(MashCommands.ParseSize("2K"), Is.EqualTo(2000));
        Assert.That(MashCommands.ParseSize("3m"), Is.EqualTo(3_000_000));
        Assert.That(MashCommands.ParseSize("1T"), Is.EqualTo(1_000_000_000_000));
        Assert.Throws<ArgumentException>(() => MashCommands.ParseSize("5X"));
        Assert.Throws<ArgumentException>(() => MashCommands.ParseSize("0"));
        Assert.Throws<ArgumentException>(() => MashCommands.ParseSize("1.5K"));
    }

    [Test]
    public void SketchFile_EqualsMash_ByteForByte()
    {
        var rows = Oracle.Value;
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(100));
        var failures = new List<string>();
        int multiSegment = 0, reads = 0, oldList = 0, newList = 0, individual = 0, list = 0, wide = 0;
        foreach (var r in rows)
        {
            byte[] got;
            try
            {
                got = Make(r, out _);
            }
            catch (ArgumentException e)
            {
                failures.Add($"row {r.Id}: {e.Message}");
                continue;
            }

            // a read-set sketch holds the multiplicities in the iteration order of Mash's robin_hood::unordered_map (A5-10b3, F194: RobinHoodMap): every file is compared byte for byte
            if (!got.AsSpan().SequenceEqual(r.Msh))
                failures.Add($"row {r.Id} {string.Join(' ', r.Args)}: expected {r.Msh.Length} bytes, got {got.Length}; first difference at {FirstDifference(got, r.Msh)}; {Describe(MashSketchFile.Read(got), MashSketchFile.Read(r.Msh))}");
            if (BitConverter.ToUInt32(r.Msh) > 0)
                multiSegment++;
            if (r.Args.Contains("-r"))
                reads++;
            if (r.Args.Contains("-i"))
                individual++;
            if (r.Args.Contains("-l"))
                list++;
            int seedIndex = r.Args.IndexOf("-S");
            if (seedIndex >= 0 && r.Args[seedIndex + 1] == "42")
                oldList++;
            else
                newList++;
            if (int.Parse(r.Args[r.Args.IndexOf("-k") + 1], CultureInfo.InvariantCulture) > 16)
                wide++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(5)));
        Assert.That(multiSegment, Is.GreaterThanOrEqualTo(15), "messages with several segments (far pointers)");
        Assert.That(reads, Is.GreaterThanOrEqualTo(20));
        Assert.That(individual, Is.GreaterThanOrEqualTo(10));
        Assert.That(list, Is.GreaterThanOrEqualTo(8));
        Assert.That(oldList, Is.GreaterThanOrEqualTo(30));
        Assert.That(newList, Is.GreaterThanOrEqualTo(30));
        Assert.That(wide, Is.GreaterThanOrEqualTo(40));
    }

    // which field of which sketch differs (a diagnosis for a failing row)
    private static string Describe(MashSketchFileData got, MashSketchFileData expected)
    {
        if (got.KmerSize != expected.KmerSize || got.HashSeed != expected.HashSeed || got.MinHashesPerWindow != expected.MinHashesPerWindow || got.Noncanonical != expected.Noncanonical)
            return "parameters differ";
        if (got.References.Count != expected.References.Count)
            return $"{got.References.Count} sketches instead of {expected.References.Count}";
        for (int i = 0; i < got.References.Count; i++)
        {
            var g = got.References[i];
            var e = expected.References[i];
            if (g.Name != e.Name || g.Comment != e.Comment || g.Length != e.Length)
                return $"sketch {i}: name / comment / length '{g.Name}' '{g.Comment}' {g.Length} vs '{e.Name}' '{e.Comment}' {e.Length}";
            if (!g.Hashes.SequenceEqual(e.Hashes))
                return $"sketch {i}: hashes differ ({g.Hashes.Count} vs {e.Hashes.Count})";
            if (!(g.Counts ?? []).SequenceEqual(e.Counts ?? []))
            {
                var gc = g.Counts ?? [];
                var ec = e.Counts ?? [];
                int d = 0;
                while (d < Math.Min(gc.Count, ec.Count) && gc[d] == ec[d]) d++;
                return $"sketch {i}: counts differ ({gc.Count} vs {ec.Count}) from index {d}: {(d < gc.Count ? gc[d] : 0)} vs {(d < ec.Count ? ec[d] : 0)}; hashes {string.Join(',', g.Hashes.Take(12))}; ours {string.Join(',', gc.Take(12))} expected {string.Join(',', ec.Take(12))}";
            }
        }

        return "same content (layout differs)";
    }

    private static string FirstDifference(byte[] a, byte[] b)
    {
        for (int i = 0; i < Math.Min(a.Length, b.Length); i++)
            if (a[i] != b[i])
                return $"byte {i} (word {i / 8})";
        return "the length";
    }

    [Test]
    public void ReadThenWrite_ReproducesTheFile()
    {
        // the reader follows far pointers; the content written again is the file
        foreach (var r in Oracle.Value)
        {
            var data = MashSketchFile.Read(r.Msh);
            Assert.That(MashSketchFile.Write(data), Is.EqualTo(r.Msh), $"row {r.Id}");
        }
    }

    [Test]
    public void Info_EqualsMash_InEveryOutputMode()
    {
        var failures = new List<string>();
        int counted = 0, errors = 0;
        foreach (var r in Oracle.Value)
        {
            var data = MashSketchFile.Read(r.Msh);
            foreach (var (key, options) in new[]
            {
                ("default", new MashCommands.MashInfoOptions()), ("t", new MashCommands.MashInfoOptions(Tabular: true)), ("H", new MashCommands.MashInfoOptions(HeaderOnly: true)),
                ("d", new MashCommands.MashInfoOptions(Dump: true)), ("c", new MashCommands.MashInfoOptions(Counts: true)),
            })
            {
                string expected = r.Info[key];
                string got;
                try
                {
                    got = MashCommands.Info(data, options);
                }
                catch (ArgumentException e)
                {
                    got = "ERROR: " + e.Message;   // the program prints "ERROR: " + the text on the error stream
                    errors++;
                }

                if (key == "c" && !got.StartsWith("ERROR", StringComparison.Ordinal))
                    counted++;
                if (got != expected)
                    failures.Add($"row {r.Id} info {key}\n--- expected\n{expected}--- actual\n{got}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count}:\n" + string.Join("\n", failures.Take(2)));
        Assert.That(counted, Is.GreaterThanOrEqualTo(15), "histograms");
        Assert.That(errors, Is.GreaterThanOrEqualTo(50), "-c on a sketch without counts is refused");
    }

    [Test]
    public void IncompatibleInfoOptions_AreRefused()
    {
        var data = MashSketchFile.Read(Oracle.Value[0].Msh);
        Assert.Throws<ArgumentException>(() => MashCommands.Info(data, new MashCommands.MashInfoOptions(HeaderOnly: true, Tabular: true)));
        Assert.Throws<ArgumentException>(() => MashCommands.Info(data, new MashCommands.MashInfoOptions(Dump: true, Counts: true)));
    }
}
