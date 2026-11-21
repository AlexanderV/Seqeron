using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-10b (B09 F185): <c>mash screen</c> of Mash 2.3 (<see cref="MashCommands.Screen"/>: the hashes of the query sketches counted in the windows of a mixture of reads or contigs; identity, shared-hashes, median-multiplicity,
/// p-value, query-ID, comment; <c>-w</c> winner-takes-all, <c>-i</c> minimum identity, <c>-v</c> maximum p-value; query sketches made with <c>-k -s -S -n -i</c>) against the installed mash 2.3 on 240 random screens
/// (TestData/Alignment/mash_screen_oracle.jsonl.gz, scripts/alignment/generate_mash_screen_oracle.py), the whole stdout compared as TEXT.
/// </summary>
[TestFixture]
public class MashScreen_Oracle_Tests
{
    private static List<MashInput> Inputs(JsonElement array) => [.. array.EnumerateArray().Select(f => new MashInput(f.GetProperty("name").GetString()!,
        [.. f.GetProperty("records").EnumerateArray().Select(r => new MashRecord(r.GetProperty("id").GetString()!, r.GetProperty("comment").GetString()!, r.GetProperty("seq").GetString()!))]))];

    [Test]
    public void Screen_EqualsMash()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.mash_screen_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        int rows = 0, withOutput = 0, winning = 0;
        var failures = new List<string>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            string expected = j.GetProperty("expected").GetString()!;
            string got = MashCommands.Screen(Inputs(j.GetProperty("queries")), Inputs(j.GetProperty("mixtures")), j.GetProperty("k").GetInt32(), j.GetProperty("s").GetInt32(), (uint)j.GetProperty("seed").GetInt64(),
                !j.GetProperty("noncanonical").GetBoolean(), j.GetProperty("individual").GetBoolean(), j.GetProperty("winning").GetBoolean(), j.GetProperty("identity").GetDouble(), j.GetProperty("pvalue").GetDouble());
            if (got != expected)
                failures.Add($"row {rows} k={j.GetProperty("k")} s={j.GetProperty("s")}\n--- expected\n{expected}--- actual\n{got}");
            rows++;
            if (expected.Length > 0)
                withOutput++;
            if (j.GetProperty("winning").GetBoolean())
                winning++;
        }

        Assert.That(rows, Is.GreaterThanOrEqualTo(200));
        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(2)));
        Assert.That(withOutput, Is.GreaterThanOrEqualTo(rows / 2), "screens with output");
        Assert.That(winning, Is.GreaterThanOrEqualTo(40), "winner-takes-all rows");
    }

    [Test]
    public void Screen_Example_QueryContainedInReads()
    {
        // mash sketch -k 8 -s 100 g.fa; mash screen g.msh m.fa (installed 2.3): every hash of the query is in the reads, once
        var q = new MashInput("g.fa", [new MashRecord("g", "", "ACGTTGCATGCCGATAGCTAGGATCCATGCATGCAAGT")]);
        var m = new MashInput("m.fa", [new MashRecord("r1", "", "ACGTTGCATGCCGATAGCTAGGATCC"), new MashRecord("r2", "", "GATAGCTAGGATCCATGCATGCAAGT")]);
        Assert.That(MashCommands.Screen([q], [m], 8, 100), Is.EqualTo("1\t30/30\t1\t6.59538e-101\tg.fa\tg \n"));
    }

    [Test]
    public void Screen_RefusesWhatTheProgramRefuses()
    {
        var q = new MashInput("q.fa", [new MashRecord("q", "", "ACGTACGTACGTAGCTAGCT")]);
        Assert.Throws<ArgumentException>(() => MashCommands.Screen([q], [new MashInput("m.fa", [])]), "no mixture record");
        Assert.Throws<ArgumentException>(() => MashCommands.Screen([new MashInput("s.fa", [new MashRecord("s", "", "ACG")])], [q]), "query shorter than k");
        Assert.Throws<ArgumentOutOfRangeException>(() => MashCommands.Screen([q], [q], 21, 1000, 42, true, false, false, 2.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => MashCommands.Screen([q], [q], 33));
    }
}
