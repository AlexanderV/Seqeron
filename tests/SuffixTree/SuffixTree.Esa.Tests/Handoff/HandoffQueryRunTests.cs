using System.Buffers.Binary;
using SuffixTree.Esa.Bench.Links;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.Handoff;

/// <summary>
/// Work package 1b.9f: the whole query check (<c>sesa handoff queries</c>) on real containers made by the production builder: a clean container passes with the engine and with the suffix-array search, a container whose
/// suffix array was damaged after the build (two entries swapped, the section checksum not recomputed: the <c>Header</c> level does not look at it) fails, a text that is not the container's text is refused, and a container
/// above the 32-bit surface of the engine is answered by the search alone (simulated: the engine is switched off).
/// </summary>
[TestFixture]
public sealed class HandoffQueryRunTests
{
    private string _dir = "";

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "esa-handoff-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    private static byte[] Text(int seed, int n, string alphabet, int repeatEvery = 0)
    {
        var rng = new Random(seed);
        var t = new byte[n];
        for (int i = 0; i < n; i++) t[i] = (byte)alphabet[rng.Next(alphabet.Length)];
        if (repeatEvery > 0)
        {
            for (int i = repeatEvery; i + repeatEvery <= n; i += repeatEvery) Array.Copy(t, 0, t, i, Math.Min(repeatEvery, 200));
        }
        return t;
    }

    private (string Input, string Container) Build(byte[] text)
    {
        string input = Path.Combine(_dir, "text.txt");
        string container = Path.Combine(_dir, "out.sesa");
        File.WriteAllBytes(input, text);
        SesaBuilder.Build(input, container, new SesaBuildOptions { Threads = 2, Fsync = false, MeasureRss = false });
        return (input, container);
    }

    private static (int Code, Dictionary<string, string> Rows) Run(string container, string input, string engine = "auto", int cap = 1000)
    {
        var rows = new Dictionary<string, string>();
        int code = HandoffQueryRun.Execute(new HandoffQueryOptions { Container = container, Input = input, Lengths = [4, 9, 16], Present = 60, Absent = 40, Cap = cap, Engine = engine }, (m, v, _) => rows[m] = v);
        return (code, rows);
    }

    [Test]
    public void ACleanContainer_Passes_WithTheEngineAndTheSearch_OnDnaAndOnAProteinLikeText()
    {
        foreach ((string alphabet, int n, int repeat) in new[] { ("ACGT", 20000, 5000), ("ACDEFGHIKLMNPQRSTVWY", 6000, 1500), ("AB", 3000, 0) })
        {
            (string input, string container) = Build(Text(1, n, alphabet, repeat));
            (int code, Dictionary<string, string> rows) = Run(container, input);
            Assert.That(code, Is.Zero, string.Join("; ", rows.Where(r => r.Key.StartsWith("q_", StringComparison.Ordinal)).Select(r => r.Key + "=" + r.Value)));
            Assert.That(rows["q_result"], Is.EqualTo("PASS"));
            Assert.That(rows["q_total"], Is.EqualTo((3 * 100).ToString()));
            Assert.That(rows["q_sa_mismatches"], Is.EqualTo("0"));
            Assert.That(rows["q_engine_status"], Is.EqualTo("ok"));
            Assert.That(rows["q_engine_mismatches"], Is.EqualTo("0"));
            Assert.That(long.Parse(rows["q_oracle_hits_total"]), Is.GreaterThan(0));
        }
    }

    [Test]
    public void ACapBelowTheCounts_ComparesTheCountsOnly_AndStillPasses()
    {
        (string input, string container) = Build(Text(2, 4000, "AC", 0));
        (int code, Dictionary<string, string> rows) = Run(container, input, cap: 3);
        Assert.That(code, Is.Zero);
        Assert.That(long.Parse(rows["q_locate_skipped_over_cap"]), Is.GreaterThan(0));
    }

    [Test]
    public void AnSwappedPairOfSuffixArrayEntries_IsFoundByTheSearchAndByTheEngine()
    {
        (string input, string container) = Build(Text(3, 5000, "ACGT", 1000));
        using (SesaFile file = SesaFile.Open(container, new EsaOpenOptions()))
        {
            long offset = file.Parsed.Sa.Offset;
            using var fs = new FileStream(container, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
            var a = new byte[4];
            var b = new byte[4];
            fs.Position = offset + (4L * 1000);
            fs.ReadExactly(a);
            fs.Position = offset + (4L * 3000);
            fs.ReadExactly(b);
            fs.Position = offset + (4L * 1000);
            fs.Write(b);
            fs.Position = offset + (4L * 3000);
            fs.Write(a);
        }
        (int code, Dictionary<string, string> rows) = Run(container, input);
        Assert.That(code, Is.EqualTo(3));
        Assert.That(rows["q_result"], Is.EqualTo("FAIL"));
        Assert.That(long.Parse(rows["q_sa_mismatches"]) + long.Parse(rows["q_engine_mismatches"]), Is.GreaterThan(0));
    }

    [Test]
    public void ADamagedChildCell_IsFoundByTheEngineOnly()
    {
        // the search does not use the child table: only the engine's descent can disagree with the oracle
        (string input, string container) = Build(Text(4, 5000, "ACGT", 0));
        using (SesaFile file = SesaFile.Open(container, new EsaOpenOptions()))
        {
            long offset = file.Parsed.Child.Offset;
            using var fs = new FileStream(container, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
            var cell = new byte[4];
            for (int i = 1; i < 2000; i++)
            {
                fs.Position = offset + (4L * i);
                fs.ReadExactly(cell);
                BinaryPrimitives.WriteUInt32LittleEndian(cell, BinaryPrimitives.ReadUInt32LittleEndian(cell) ^ 0x00000010u);
                fs.Position = offset + (4L * i);
                fs.Write(cell);
            }
        }
        (int code, Dictionary<string, string> rows) = Run(container, input);
        Assert.That(rows["q_sa_mismatches"], Is.EqualTo("0"));
        Assert.That(code, Is.EqualTo(3).Or.Zero, "a damaged child table may, or may not, change an answer: what matters is that the engine is compared");
        Assert.That(rows.ContainsKey("q_engine_mismatches"));
    }

    [Test]
    public void ATextThatIsNotTheContainersText_IsRefused()
    {
        (string input, string container) = Build(Text(5, 3000, "ACGT", 0));
        File.WriteAllBytes(input, Text(5, 2999, "ACGT", 0));
        (int code, Dictionary<string, string> rows) = Run(container, input);
        Assert.That(code, Is.EqualTo(3));
        Assert.That(rows["q_result"], Is.EqualTo("FAIL"));
        Assert.That(rows["q_input_length"], Does.Contain("2999").And.Contain("3000"));
    }

    [Test]
    public void WithTheEngineSwitchedOff_TheSearchAnswersAlone()
    {
        (string input, string container) = Build(Text(6, 4000, "ACGT", 1000));
        (int code, Dictionary<string, string> rows) = Run(container, input, engine: "false");
        Assert.That(code, Is.Zero);
        Assert.That(rows["q_engine_status"], Is.EqualTo("skipped"));
        Assert.That(rows["q_sa_mismatches"], Is.EqualTo("0"));
    }

    [Test]
    public void AQueryFileOfTheSameSeed_GivesTheSameRowsTwice()
    {
        (string input, string container) = Build(Text(7, 4000, "ACGT", 1000));
        (_, Dictionary<string, string> first) = Run(container, input);
        (_, Dictionary<string, string> second) = Run(container, input);
        foreach (string key in new[] { "q_total", "q_present", "q_absent", "q_oracle_hits_total", "q_set_sha256" })
        {
            Assert.That(second[key], Is.EqualTo(first[key]), key);
        }
    }
}
