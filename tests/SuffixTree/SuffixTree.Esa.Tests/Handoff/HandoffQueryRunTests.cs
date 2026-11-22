using System.Buffers.Binary;
using System.Reflection;
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
    public void AReversedBlockOfSuffixArrayEntries_IsFoundByTheSearchAndByTheEngine()
    {
        // one damaged entry is the business of e8-verify (the order check); a query set finds damage that a query can reach, so the damage here is a block of 1,000 of the 5,001 entries
        (string input, string container) = Build(Text(3, 5000, "ACGT", 1000));
        using (SesaFile file = SesaFile.Open(container, new EsaOpenOptions()))
        {
            long offset = file.Parsed.Sa.Offset + (4L * 1000);
            using var fs = new FileStream(container, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
            var block = new byte[4 * 1000];
            fs.Position = offset;
            fs.ReadExactly(block);
            var reversed = new byte[block.Length];
            for (int i = 0; i < 1000; i++) Array.Copy(block, 4 * (999 - i), reversed, 4 * i, 4);
            fs.Position = offset;
            fs.Write(reversed);
        }
        (int code, Dictionary<string, string> rows) = Run(container, input);
        Assert.That(code, Is.EqualTo(3));
        Assert.That(rows["q_result"], Is.EqualTo("FAIL"));
        Assert.That(long.Parse(rows["q_sa_mismatches"]), Is.GreaterThan(0), "the search");
        Assert.That(long.Parse(rows["q_engine_mismatches"]), Is.GreaterThan(0), "the engine reads the positions of a locate from the same array");
    }

    [Test]
    public void ADamagedChildTable_IsFoundByTheEngineAndNotByTheSearch()
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
        Assert.That(code, Is.EqualTo(3));
        Assert.That(long.Parse(rows["q_engine_mismatches"]), Is.GreaterThan(0));
        Assert.That(long.Parse(rows["q_engine_refused"]), Is.GreaterThan(0), "the Checked mode refuses a damaged child table instead of answering");
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
        Assert.That(rows.ContainsKey("q_total"), Is.False, "refused before any query is cut (the audit's mutant that lets the run go on was not seen: its result is FAIL as well)");
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

    private static void ReverseSaBlock(string container, int from, int count)
    {
        using SesaFile file = SesaFile.Open(container, new EsaOpenOptions());
        long offset = file.Parsed.Sa.Offset + (4L * from);
        using var fs = new FileStream(container, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
        var block = new byte[4 * count];
        fs.Position = offset;
        fs.ReadExactly(block);
        var reversed = new byte[block.Length];
        for (int i = 0; i < count; i++) Array.Copy(block, 4 * (count - 1 - i), reversed, 4 * i, 4);
        fs.Position = offset;
        fs.Write(reversed);
    }

    [Test]
    public void ADamagedSuffixArray_FailsTheRun_ByTheSearchAlone_WhenTheEngineIsSwitchedOff()
    {
        // the audit's mutant "the mismatches of the search do not fail the run" was hidden by the engine, which reads the same array and fails the run by itself
        (string input, string container) = Build(Text(3, 5000, "ACGT", 1000));
        ReverseSaBlock(container, 1000, 1000);
        (int code, Dictionary<string, string> rows) = Run(container, input, engine: "false");
        Assert.That(long.Parse(rows["q_sa_mismatches"]), Is.GreaterThan(0));
        Assert.That(rows["q_result"], Is.EqualTo("FAIL"));
        Assert.That(code, Is.EqualTo(3));
    }

    [Test]
    public void ExactlyOneWrongAnswer_FailsTheRun()
    {
        // one wrong entry of the suffix array in the range of a query that occurs twice: the count and the flag stay right, the locate set has the position of the other occurrence twice: ONE mismatch
        // (the audit's mutant "one mismatch is tolerated" passed such a run)
        byte[] text = Text(8, 4000, "ACGT", 0);
        Array.Copy(text, 0, text, 2000, 200);
        string textPath = Path.Combine(_dir, "probe.txt");
        File.WriteAllBytes(textPath, text);
        ulong seed = 0;
        long a = -1, b = -1;
        for (ulong s = 1; s < 500 && a < 0; s++)
        {
            IReadOnlyList<HandoffQuery> set = HandoffQuerySet.Build(text.Length, (pos, len) => text.AsSpan((int)pos, len).ToArray(), "ACGT"u8, [16], 1, 0, s);
            ScanResult r = RollingScan.Run(new MemoryStream(text), set, 10);
            if (r.Counts[0] == 2) (seed, a, b) = (s, r.Positions[0][0], r.Positions[0][1]);
        }
        Assert.That(a, Is.GreaterThanOrEqualTo(0), "a seed whose one query occurs twice");
        (string input, string container) = Build(text);
        using (SesaFile file = SesaFile.Open(container, new EsaOpenOptions()))
        {
            using var fs = new FileStream(container, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
            var cell = new byte[4];
            for (long rank = 0; rank <= text.Length; rank++)
            {
                fs.Position = file.Parsed.Sa.Offset + (4 * rank);
                fs.ReadExactly(cell);
                if (BinaryPrimitives.ReadUInt32LittleEndian(cell) != (uint)a) continue;
                BinaryPrimitives.WriteUInt32LittleEndian(cell, (uint)b);
                fs.Position = file.Parsed.Sa.Offset + (4 * rank);
                fs.Write(cell);
                break;
            }
        }
        var rows = new Dictionary<string, string>();
        int code = HandoffQueryRun.Execute(new HandoffQueryOptions { Container = container, Input = input, Lengths = [16], Present = 1, Absent = 0, Seed = seed, Engine = "false" }, (m, v, _) => rows[m] = v);
        Assert.That(rows["q_sa_mismatches"], Is.EqualTo("1"), string.Join("; ", rows.Where(r => r.Key.StartsWith("q_", StringComparison.Ordinal)).Select(r => r.Key + "=" + r.Value)));
        Assert.That(rows["q_result"], Is.EqualTo("FAIL"));
        Assert.That(code, Is.EqualTo(3));
    }

    /// <summary>A tree that answers like the real one except that <c>Contains</c> says the opposite.</summary>
    public class ContainsLies : DispatchProxy
    {
        public ISuffixTree? Inner { get; set; }

        protected override object? Invoke(MethodInfo? targetMethod, object?[]? args)
        {
            try
            {
                object? r = targetMethod!.Invoke(Inner, args);
                return targetMethod.Name == nameof(ISuffixTree.Contains) && r is bool b ? !b : r;
            }
            catch (TargetInvocationException e) when (e.InnerException is not null)
            {
                throw e.InnerException;
            }
        }
    }

    [Test]
    public void AnEngineWhoseContainsIsWrong_ButWhoseCountIsRight_FailsTheRun()
    {
        // the engine's Contains is asked on its own (the audit's mutant that took "exists" from the count passed this engine)
        (string input, string container) = Build(Text(9, 4000, "ACGT", 1000));
        ISuffixTree inner = SuffixIndexFactory.Open(container, new EsaOpenOptions());
        try
        {
            var rows = new Dictionary<string, string>();
            int code = HandoffQueryRun.Execute(new HandoffQueryOptions
            {
                Container = container,
                Input = input,
                Lengths = [4, 9, 16],
                Present = 60,
                Absent = 40,
                OpenEngine = _ =>
                {
                    ISuffixTree liar = DispatchProxy.Create<ISuffixTree, ContainsLies>();
                    ((ContainsLies)liar).Inner = inner;
                    return liar;
                },
            }, (m, v, _) => rows[m] = v);
            Assert.That(rows["q_sa_mismatches"], Is.EqualTo("0"));
            Assert.That(long.Parse(rows["q_engine_mismatches"]), Is.EqualTo(300), "every query: the flag is the opposite of the truth");
            Assert.That(code, Is.EqualTo(3));
        }
        finally
        {
            (inner as IDisposable)?.Dispose();
        }
    }
}
