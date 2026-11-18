using System.Security.Cryptography;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The stage of work package 1b.7c as a whole (PH4 base mode fused with PH5 through the B3 hook, on a given suffix-array file): the contract against the in-RAM oracles on small texts and on the families, every
/// width, the chunk / block / window boundaries, the thread counts, both codecs, the failure modes, the counters and the platform without <c>/proc</c>.
/// </summary>
[TestFixture]
public class LcpChildBuilderTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private string Scratch => _dir.Sub("scratch");

    private LcpChildOptions Tiny(int v, int threads = 1, ChildCodec codec = ChildCodec.Packed, int? sigma = null) => LcpChildRun.Options(Scratch, LcpChildRun.TinyTuning(v), threads, codec, sigma);

    [Test]
    public void SmallTexts_AllOutputs_EqualTheInRamStructures_BothCodecs()
    {
        int v = 0;
        foreach (string text in new[] { "a", "ab", "ba", "aa", "banana", "mississippi", "abracadabra", "aaaaaaaaaaaaaaaaaa", "abababababababab", "abcabcabcabcabc", "aabaabaabaab", StringFamilies.Fibonacci(80), StringFamilies.ThueMorse(70), StringFamilies.RandomDna(100, 4), StringFamilies.RandomProtein(90, 5) })
        {
            LcpChildCase c = LcpChildCase.Of(text);
            foreach (ChildCodec codec in new[] { ChildCodec.Packed, ChildCodec.Tagged })
            {
                using var dir = new TempDir();
                LcpChildRun.Outcome o = LcpChildRun.Run(dir, c, LcpChildRun.Options(dir.Sub("scratch"), LcpChildRun.TinyTuning(v), 1 + (v % 3), codec, c.Sigma), stats: true);
                LcpChildRun.Compare(c, o, codec, $"'{(text.Length > 30 ? text[..30] : text)}' {codec} variant {v}");
                v++;
            }
        }
    }

    [Test]
    public void TheEmptyText_And_TextsOfOneAndTwoSymbols_Work()
    {
        foreach (string text in new[] { string.Empty, "x", "xx", "xy", "yx" })
        {
            LcpChildCase c = LcpChildCase.Of(text);
            LcpChildRun.Outcome o = LcpChildRun.Run(_dir, c, Tiny(1, sigma: c.Sigma), tag: text.Length.ToString() + text);
            LcpChildRun.Compare(c, o, ChildCodec.Packed, $"'{text}'");
            Assert.That(new FileInfo(o.Files.Child).Length, Is.EqualTo(4 * (text.Length + 1)));
            Assert.That(o.Result.Scan.NodeCount, Is.EqualTo(c.Statistics.NodeCount));
        }
    }

    [Test]
    public void Families_AtEveryTinyVariant_EqualTheInRamStructures()
    {
        int v = 0;
        foreach (int n in new[] { 600, 3000 })
        {
            foreach ((string name, string text) in LcpChildRun.Families(n))
            {
                LcpChildCase c = LcpChildCase.Of(text);
                ChildCodec codec = v % 2 == 0 ? ChildCodec.Packed : ChildCodec.Tagged;
                using var dir = new TempDir();
                LcpChildTuning t = LcpChildRun.TinyTuning(v) with { LcpChunk = 17 + (v * 29 % 400) };
                LcpChildRun.Outcome o = LcpChildRun.Run(dir, c, LcpChildRun.Options(dir.Sub("scratch"), t, 1 + (v % 4), codec, c.Sigma));
                LcpChildRun.Compare(c, o, codec, $"{name} {n} variant {v}");
                v++;
            }
        }
    }

    [Test]
    public void Families_WithTheDefaultParameters_EqualTheInRamStructures()
    {
        foreach ((string name, string text) in LcpChildRun.Families(60_000))
        {
            LcpChildCase c = LcpChildCase.Of(text);
            using var dir = new TempDir();
            LcpChildRun.Outcome o = LcpChildRun.Run(dir, c, new LcpChildOptions { ScratchDir = dir.Sub("scratch"), Sigma = c.Sigma });
            LcpChildRun.Compare(c, o, ChildCodec.Packed, $"{name} 60000, default parameters");
            Assert.That(o.Result.Parameters.LcpChunk, Is.EqualTo(8_000_000));
        }
    }

    [Test]
    public void EveryWidth_IsMet_WithItsOverflow_AndTheChildTableOfTheSameText()
    {
        var kinds = new HashSet<LcpKind>();
        string u8Overflow = StringFamilies.RandomDna(20_000, 5) + StringFamilies.RandomDna(600, 11) + StringFamilies.RandomDna(100, 6) + StringFamilies.RandomDna(600, 11);
        string u16 = StringFamilies.RandomDna(20_000, 7) + StringFamilies.RandomDna(70_000, 8) + StringFamilies.RandomDna(333, 9) + StringFamilies.RandomDna(70_000, 8);
        foreach ((string name, string text) in new[] { ("u8", StringFamilies.RandomDna(20_000, 1)), ("u8 + overflow", u8Overflow), ("u16 + overflow", u16), ("u32", new string('a', 4000)) })
        {
            LcpChildCase c = LcpChildCase.Of(text);
            kinds.Add(c.Lcp.Kind);
            using var dir = new TempDir();
            LcpChildRun.Outcome o = LcpChildRun.Run(dir, c, LcpChildRun.Options(dir.Sub("scratch"), new LcpChildTuning { LcpChunk = 5000, BlockEntries = 256, ChildWindowBytes = 4096 }, 2, ChildCodec.Packed, c.Sigma));
            LcpChildRun.Compare(c, o, ChildCodec.Packed, name);
            if (name == "u8 + overflow") Assert.That(o.Result.Lcp.OverflowCount, Is.GreaterThan(0).And.LessThan(2000));
            if (name == "u16 + overflow") Assert.That((o.Result.Lcp.Kind, o.Result.Lcp.OverflowCount), Is.EqualTo((LcpKind.U16, (long)c.Lcp.OverflowCount)));
            if (name == "u32") Assert.That(new FileInfo(o.Files.Overflow).Length, Is.Zero);
        }
        Assert.That(kinds, Is.EquivalentTo(new[] { LcpKind.U8, LcpKind.U16, LcpKind.U32 }));
    }

    [Test]
    public void ChunkBoundaries_EveryNAroundR_AndAround2R()
    {
        int v = 0;
        foreach (int r in new[] { 3, 4, 7, 16 })
        {
            foreach (int n in new[] { r - 1, r, r + 1, (2 * r) - 1, 2 * r, (2 * r) + 1, (5 * r) - 1, 5 * r, (5 * r) + 1 })
            {
                if (n < 1) continue;
                foreach (string text in new[] { StringFamilies.RandomDna(n, (ulong)(r * 100) + (ulong)n), new string('a', n), StringFamilies.Fibonacci(n) })
                {
                    LcpChildCase c = LcpChildCase.Of(text);
                    using var dir = new TempDir();
                    LcpChildTuning t = LcpChildRun.TinyTuning(v) with { LcpChunk = r };
                    LcpChildRun.Outcome o = LcpChildRun.Run(dir, c, LcpChildRun.Options(dir.Sub("scratch"), t, 1 + (v % 3), ChildCodec.Packed, c.Sigma));
                    LcpChildRun.Compare(c, o, ChildCodec.Packed, $"R={r} n={n} variant {v}");
                    Assert.That(o.Result.Lcp.TextChunks, Is.EqualTo((n + r - 1) / r));
                    Assert.That(o.Result.Lcp.RankChunks, Is.EqualTo((n + r) / r));
                    v++;
                }
            }
        }
    }

    [Test]
    public void SpillBoundaries_ADepthAroundTwoBlocks_SpillsExactlyWhenItShould()
    {
        foreach (int block in new[] { 2, 3, 6, 17 })
        {
            foreach (int n in new[] { (2 * block) - 1, 2 * block, (2 * block) + 1, (2 * block) + 2, 3 * block, (3 * block) + 1, 9 * block })
            {
                LcpChildCase c = LcpChildCase.Of(new string('a', n));
                using var dir = new TempDir();
                LcpChildTuning t = new() { LcpChunk = 10, BlockEntries = block, ChildWindowBytes = 32, ChildSlots = 2, SpillFlushBytes = 1 };
                LcpChildRun.Outcome o = LcpChildRun.Run(dir, c, LcpChildRun.Options(dir.Sub("scratch"), t, 2, ChildCodec.Packed, c.Sigma));
                LcpChildRun.Compare(c, o, ChildCodec.Packed, $"a^{n} block {block}");
                Assert.That(o.Result.Scan.StackPeak, Is.EqualTo(n));
                Assert.That(o.Result.Scan.Spill.BlocksWritten > 0, Is.EqualTo(n > 2 * block), $"a^{n} block {block}");
            }
        }
    }

    [Test]
    public void TheOutput_IsIdentical_For1_2_3_And4Threads_AndForEveryWorkerCount()
    {
        foreach ((string name, string text) in new[] { ("random", StringFamilies.RandomDna(2500, 3)), ("a^n", new string('a', 1500)), ("fibonacci", StringFamilies.Fibonacci(2000)), ("repeat-rich", LcpChildRun.RepeatRich(2500, 9)) })
        {
            LcpChildCase c = LcpChildCase.Of(text);
            string[]? reference = null;
            int v = 0;
            foreach (int threads in new[] { 1, 2, 3, 4 })
            {
                foreach (int workers in new[] { 1, 2, 4 })
                {
                    using var dir = new TempDir();
                    LcpChildTuning t = LcpChildRun.TinyTuning(v++) with { LcpChunk = 60, LcpWorkers = workers };
                    LcpChildRun.Outcome o = LcpChildRun.Run(dir, c, LcpChildRun.Options(dir.Sub("scratch"), t, threads, ChildCodec.Packed, c.Sigma));
                    LcpChildRun.Compare(c, o, ChildCodec.Packed, $"{name} threads {threads} workers {workers}");
                    string[] h = LcpChildRun.Hashes(o);
                    reference ??= h;
                    Assert.That(h, Is.EqualTo(reference), $"{name} threads {threads} workers {workers}");
                }
            }
        }
    }

    [Test]
    public void TheDefaultRule_Runs_TwoWorkers_OfFour()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(40_000, 4));
        using var dir = new TempDir();
        LcpChildRun.Outcome o = LcpChildRun.Run(dir, c, new LcpChildOptions { ScratchDir = dir.Sub("scratch"), Threads = 4, Sigma = c.Sigma, Tuning = new LcpChildTuning { LcpChunk = 1000 } });
        LcpChildRun.Compare(c, o, ChildCodec.Packed, "rule workers");
        Assert.That(o.Result.Parameters.LcpWorkers, Is.EqualTo(2));
        Assert.That(o.Result.Lcp.Workers, Is.EqualTo(2));
        Assert.That(o.Result.Parameters.LcpChunk, Is.EqualTo(1000));
    }

    [Test]
    public void TheTaggedCodec_IsRefusedAbove2Pow30Minus1_AndAnyTextAboveTheDesignLimit_BeforeAnythingIsTouched()
    {
        string text = _dir.File("big.pk");
        long n = 1L << 30;
        using (var f = new FileStream(text, FileMode.Create)) f.SetLength(PackedTextFormat.TotalBytes(n, 1)); // sparse: 2^30 symbols of one bit
        using MappedTextReader reader = MappedTextReader.Open(text, n, 1);
        string outDir = _dir.Sub("out");
        Directory.CreateDirectory(outDir);
        var files = new LcpChildOutputs(Path.Combine(outDir, "a"), Path.Combine(outDir, "b"), Path.Combine(outDir, "c"));
        string sa = _dir.File("missing.sa"); // never opened: the codec check comes first
        Assert.Throws<NotSupportedException>(() => LcpChildBuilder.Build(sa, reader, files, new LcpChildOptions { ScratchDir = Scratch, Codec = ChildCodec.Tagged }));
        Assert.That(Directory.GetFileSystemEntries(outDir), Is.Empty);
        Assert.That(Directory.Exists(Scratch), Is.False);

        string huge = _dir.File("huge.pk");
        long nh = PackedTextFormat.MaxSymbols + 1;
        using (var f = new FileStream(huge, FileMode.Create)) f.SetLength(PackedTextFormat.TotalBytes(nh, 1));
        using MappedTextReader hr = MappedTextReader.Open(huge, nh, 1);
        Assert.Throws<NotSupportedException>(() => LcpChildBuilder.Build(sa, hr, files, new LcpChildOptions { ScratchDir = Scratch }));
        Assert.That(Directory.GetFileSystemEntries(outDir), Is.Empty);
        Assert.That(Directory.Exists(Scratch), Is.False);
    }

    // ---- corrupt input and failing environment

    [Test]
    [CancelAfter(120_000)]
    public void ACorruptSuffixArray_EndsWithAnException_AndLeavesNothing()
    {
        LcpChildCase c = LcpChildCase.Of("mississippi");
        string textPath = TextOracle.WriteFile(_dir, "c.pk", c.Symbols, c.Bits);
        uint[] good = c.Sa.Select(x => (uint)x).ToArray();
        var cases = new List<(string Name, uint[] Sa)>
        {
            ("SA[0] != n", [.. good.Skip(1), 0u]),
            ("repeated entry", [.. good.Take(4), good[3], .. good.Skip(5)]),
            ("entry >= n", [.. good.Take(5), (uint)c.N, .. good.Skip(6)]),
            ("short", good[..^1]),
            ("long", [.. good, 0u]),
            ("empty", []),
        };
        int tag = 0;
        foreach ((string name, uint[] sa) in cases)
        {
            string saPath = _dir.File($"c{tag}.sa");
            FileBits.WriteUInts(saPath, sa);
            using var dir = new TempDir();
            foreach (long r in new long[] { 2, 1000 })
            {
                LcpChildTuning t = LcpChildRun.TinyTuning(tag++) with { LcpChunk = r };
                LcpChildRun.AssertFailsCleanly(dir, saPath, textPath, c.N, c.Bits, LcpChildRun.Options(dir.Sub("scratch"), t, 2, ChildCodec.Packed, c.Sigma), typeof(InvalidDataException), $"{name} R={r}");
            }
        }
    }

    [Test]
    public void ASpillThatFails_IsAnIoException_TheBuildEnds_AndEveryOutputIsRemoved()
    {
        LcpChildCase c = LcpChildCase.Of(new string('a', 200));
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        var t = new LcpChildTuning { LcpChunk = 50, BlockEntries = 2, ChildWindowBytes = 32, ChildSlots = 2 };
        string scratch = Scratch;
        LcpChildOptions o = LcpChildRun.Options(scratch, t, 2, ChildCodec.Packed, c.Sigma) with
        {
            Trace = step =>
            {
                if (step == "b3") Directory.Delete(Path.Combine(scratch, "ph5"), true); // the stack can no longer be spilled
            },
        };
        LcpChildRun.AssertFailsCleanly(_dir, sa, textPath, c.N, c.Bits, o, typeof(IOException), "spill failure");
    }

    [Test]
    public void AFullDisk_AtTheCreationOfAnyOutput_IsAnIoException_AndNoOutputIsLeft()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(300, 5));
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        for (int succeed = 0; succeed <= 2; succeed++)
        {
            var source = new FailingAfter(succeed);
            LcpChildOptions o = LcpChildRun.Options(Scratch, LcpChildRun.TinyTuning(succeed), 2, ChildCodec.Packed, c.Sigma) with { Files = new MappedFileOptions { Source = source } };
            using var dir = new TempDir();
            LcpChildRun.AssertFailsCleanly(dir, sa, textPath, c.N, c.Bits, o, typeof(IOException), $"the creation number {succeed + 1} fails");
            Assert.That(source.Calls, Is.EqualTo(succeed + 1));
        }
    }

    [Test]
    public void AFailureAfterTheOutputsAreComplete_AStatsFileThatCannotBeWritten_RemovesThemAll()
    {
        LcpChildCase c = LcpChildCase.Of("banana");
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        string outDir = _dir.Sub("out");
        Directory.CreateDirectory(outDir);
        Directory.CreateDirectory(Path.Combine(outDir, "dir.stats")); // a directory where the stats file should go
        var files = new LcpChildOutputs(Path.Combine(outDir, "o.lcp"), Path.Combine(outDir, "o.ovf"), Path.Combine(outDir, "o.cld"), Path.Combine(outDir, "dir.stats"));
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        Assert.Throws<UnauthorizedAccessException>(() => LcpChildBuilder.Build(sa, text, files, Tiny(2, sigma: c.Sigma)));
        Assert.That(Directory.GetFileSystemEntries(outDir).Select(Path.GetFileName), Is.EqualTo(new[] { "dir.stats" }), "the other three outputs are gone, the foreign directory stays");
        Assert.That(Directory.Exists(Scratch), Is.False);
    }

    [Test]
    public void ARejectedSuffixArray_TouchesNoExistingOutput_AndLeavesAForeignScratchFile()
    {
        LcpChildCase c = LcpChildCase.Of("mississippi");
        string textPath = TextOracle.WriteFile(_dir, "s.pk", c.Symbols, c.Bits);
        string sa = _dir.File("s.sa");
        FileBits.WriteUInts(sa, [1u, 2u]);
        string outDir = _dir.Sub("out");
        Directory.CreateDirectory(outDir);
        var files = new LcpChildOutputs(Path.Combine(outDir, "o.lcp"), Path.Combine(outDir, "o.ovf"), Path.Combine(outDir, "o.cld"));
        foreach (string p in new[] { files.Lcp, files.Overflow, files.Child }) File.WriteAllText(p, "stale");
        Directory.CreateDirectory(Scratch);
        File.WriteAllText(Path.Combine(Scratch, "foreign.txt"), "mine");
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        Assert.Throws<InvalidDataException>(() => LcpChildBuilder.Build(sa, text, files, Tiny(2)));
        foreach (string p in new[] { files.Lcp, files.Overflow, files.Child }) Assert.That(File.ReadAllText(p), Is.EqualTo("stale"), "no output was created yet");
        Assert.That(File.ReadAllText(Path.Combine(Scratch, "foreign.txt")), Is.EqualTo("mine"), "the scratch directory was not ours");
        Assert.That(Directory.GetFileSystemEntries(Scratch), Has.Length.EqualTo(1), "our subdirectories are gone");
    }

    [Test]
    public void StaleOutputs_AreReplaced_AndASecondRunInTheSameDirectoriesGivesTheSameFiles()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(500, 4));
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        string outDir = _dir.Sub("out");
        Directory.CreateDirectory(outDir);
        var files = new LcpChildOutputs(Path.Combine(outDir, "o.lcp"), Path.Combine(outDir, "o.ovf"), Path.Combine(outDir, "o.cld"), Path.Combine(outDir, "o.stats"));
        File.WriteAllBytes(files.Lcp, new byte[100_000]);
        File.WriteAllText(files.Stats!, "stale");
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        string[][] hashes = new string[2][];
        for (int i = 0; i < 2; i++)
        {
            LcpChildBuilder.Build(sa, text, files, Tiny(i, sigma: c.Sigma));
            hashes[i] = [.. new[] { files.Lcp, files.Overflow, files.Child, files.Stats }.Select(p => Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(p!))))];
        }
        Assert.That(hashes[1], Is.EqualTo(hashes[0]));
        Assert.That(new FileInfo(files.Lcp).Length, Is.EqualTo(c.N + 1), "the stale 100,000 bytes are gone");
    }

    [Test]
    public void TheScratchDirectory_IsCreatedByTheBuild_AndRemovedWithItsSubdirectories()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(200, 4));
        var seen = new List<string>();
        LcpChildOptions o = Tiny(3, sigma: c.Sigma) with { Trace = s => seen.Add($"{s}:{Directory.Exists(Scratch)}:{Directory.Exists(Path.Combine(Scratch, "ph4"))}:{Directory.Exists(Path.Combine(Scratch, "ph5"))}") };
        LcpChildRun.Outcome r = LcpChildRun.Run(_dir, c, o);
        Assert.That(seen, Is.EqualTo(new[] { "b1:True:True:True", "b2:True:True:True", "b3:True:True:True" }));
        LcpChildRun.Compare(c, r, ChildCodec.Packed, "scratch");
    }

    // ---- the options

    [Test]
    public void TheOptions_AreChecked()
    {
        LcpChildCase c = LcpChildCase.Of("banana");
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        string outDir = _dir.Sub("out");
        Directory.CreateDirectory(outDir);
        var files = new LcpChildOutputs(Path.Combine(outDir, "o.lcp"), Path.Combine(outDir, "o.ovf"), Path.Combine(outDir, "o.cld"));
        var withStats = files with { Stats = Path.Combine(outDir, "o.stats") };
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        Assert.Throws<ArgumentNullException>(() => LcpChildBuilder.Build(sa, null!, files, Tiny(1)));
        Assert.Throws<ArgumentNullException>(() => LcpChildBuilder.Build(sa, text, null!, Tiny(1)));
        Assert.Throws<ArgumentNullException>(() => LcpChildBuilder.Build(sa, text, files, null!));
        Assert.Throws<ArgumentException>(() => LcpChildBuilder.Build(sa, text, withStats, Tiny(1)), "the statistics file records sigma");
        Assert.Throws<ArgumentOutOfRangeException>(() => LcpChildBuilder.Build(sa, text, files, Tiny(1) with { Threads = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => LcpChildBuilder.Build(sa, text, files, Tiny(1) with { ArenaCapBytes = -1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => LcpChildBuilder.Build(sa, text, files, Tiny(1) with { ChildArenaCapBytes = -1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => LcpChildBuilder.Build(sa, text, files, Tiny(1) with { Sigma = -1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => LcpChildBuilder.Build(sa, text, files, Tiny(1) with { Tuning = new LcpChildTuning { ChildWindowBytes = 24 } }), "a window is a power of two");
        Assert.That(Directory.GetFileSystemEntries(outDir), Is.Empty);
        Assert.That(Directory.Exists(Scratch), Is.False);
    }

    // ---- statistics

    [Test]
    public void TheStatsText_HasTheKeysAndOrderOfThePrototype_AndTheFileHoldsIt()
    {
        LcpChildCase c = LcpChildCase.Of("mississippi");
        LcpChildRun.Outcome o = LcpChildRun.Run(_dir, c, Tiny(2, sigma: c.Sigma));
        string text = File.ReadAllText(o.Files.Stats!);
        Assert.That(text, Is.EqualTo(o.Result.StatsText(c.Sigma)));
        Assert.That(text.Split('\n', StringSplitOptions.RemoveEmptyEntries).Select(l => l.Split('=')[0]), Is.EqualTo(new[] { "nodecount", "internal", "lrsdepth", "lrsrank", "n", "sigma", "lcpwidth", "overflow", "maxlcp" }));
        Assert.That(text, Does.EndWith("\n"));
        Assert.That(o.Result.StatsText(7), Does.Contain("\nsigma=7\n"));
        Assert.That(text, Is.EqualTo(c.ExpectedStats(1)));
    }

    // ---- counters and the platform

    [Test]
    public void TheCounters_OfTheBuild_AreTheTupleVolume_ThePhysicalOutputs_AndTheSpill()
    {
        foreach (string text in new[] { StringFamilies.RandomDna(1500, 3), new string('b', 300), "banana" })
        {
            LcpChildCase c = LcpChildCase.Of(text);
            var counters = new IoCounters();
            using var dir = new TempDir();
            LcpChildOptions o = LcpChildRun.Options(dir.Sub("scratch"), LcpChildRun.TinyTuning(2) with { BlockEntries = 4, LcpChunk = 64 }, 2, ChildCodec.Packed, c.Sigma) with { Counters = counters };
            LcpChildRun.Outcome r = LcpChildRun.Run(dir, c, o);
            LcpChildRun.Compare(c, r, ChildCodec.Packed, text);
            long n = c.N;
            IoSnapshot io = r.Result.Io;
            Assert.That(io.ReadBytes, Is.EqualTo((4 * (n + 1)) + (12 * n) + (8 * n) + r.Result.Scan.Spill.BytesRead), "SA, Phi tuples, rank tuples, spill");
            Assert.That(io.WriteBytes, Is.EqualTo((12 * n) + (8 * n) + r.Result.Scan.Spill.BytesWritten));
            Assert.That(io.MappedWriteBytes, Is.EqualTo(new FileInfo(r.Files.Lcp).Length + new FileInfo(r.Files.Overflow).Length + (4 * r.Result.Scan.ChildWindows.Writes)));
            Assert.That(io.ScratchLiveBytes, Is.Zero);
            Assert.That(io.ScratchPeakBytes, Is.GreaterThanOrEqualTo(12 * n));
            Assert.That(counters.ReadBytes, Is.EqualTo(io.ReadBytes), "the caller's instance is the one that counted");
            Assert.That(r.Result.Scan.ChildWindows.Writes, Is.GreaterThan(0));
        }
    }

    [Test]
    public void ThePrivateCounters_AreUsed_WhenTheCallerGivesNone()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(400, 3));
        LcpChildRun.Outcome r = LcpChildRun.Run(_dir, c, Tiny(1, sigma: c.Sigma));
        Assert.That(r.Result.Io.ReadBytes, Is.GreaterThan(0));
        Assert.That(r.Result.Io.WriteBytes, Is.GreaterThan(0));
    }

    [Test]
    public void TheSpillOfAn_AToTheN_IsSixteenBytesPerEntryBeyondTwoBlocks()
    {
        const int n = 5000, block = 100;
        LcpChildCase c = LcpChildCase.Of(new string('a', n));
        using var dir = new TempDir();
        var t = new LcpChildTuning { LcpChunk = 1000, BlockEntries = block, ChildWindowBytes = 1024, ChildSlots = 4 };
        LcpChildRun.Outcome r = LcpChildRun.Run(dir, c, LcpChildRun.Options(dir.Sub("scratch"), t, 2, ChildCodec.Packed, c.Sigma));
        LcpChildRun.Compare(c, r, ChildCodec.Packed, "a^5000");
        SpillStats s = r.Result.Scan.Spill;
        Assert.That(s.BytesWritten, Is.EqualTo(((n / block) - 2) * block * 16L), "a block for every block of the chain beyond the two in RAM");
        Assert.That(s.PeakFileBytes, Is.EqualTo(s.BytesWritten));
        Assert.That(s.BytesRead, Is.EqualTo(s.BytesWritten), "the descent reads every spilled block back once");
    }

    [Test]
    public void WithoutProc_TheBuildStillWorks_AndSaysTheRssFigureIsUnavailable()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(800, 3));
        var source = new FakeProcSource(); // every file is missing, as on macOS or Windows
        LcpChildRun.Outcome r = LcpChildRun.Run(_dir, c, Tiny(2, sigma: c.Sigma) with { ProcSource = source });
        LcpChildRun.Compare(c, r, ChildCodec.Packed, "no /proc");
        Assert.That(r.Result.RssAnonPeakMb, Is.Null);
        Assert.That(r.Result.RssVerdict, Is.EqualTo(RssVerdict.Unavailable));
        Assert.That(source.Calls, Is.GreaterThan(0), "the source was asked");
    }

    [Test]
    public void AHostileProcSource_IsUnavailable_NotAnException()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(300, 3));
        var source = new FakeProcSource { Throws = new IOException("no such device") };
        LcpChildRun.Outcome r = LcpChildRun.Run(_dir, c, Tiny(2, sigma: c.Sigma) with { ProcSource = source });
        LcpChildRun.Compare(c, r, ChildCodec.Packed, "hostile /proc");
        Assert.That(r.Result.RssVerdict, Is.EqualTo(RssVerdict.Unavailable));
    }

    [Test]
    public void TheRssVerdict_IsWithin_OrExceeded_AgainstTheCap_WhenTheFigureIsThere()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(300, 3));
        foreach ((long anonKb, double cap, RssVerdict expected, double mb) in new[] { (100_000L, 320.0, RssVerdict.Within, 102.4), (312_500L, 320.0, RssVerdict.Within, 320.0), (312_501L, 320.0, RssVerdict.Exceeded, 320.00102), (400_000L, 320.0, RssVerdict.Exceeded, 409.6), (400_000L, 500.0, RssVerdict.Within, 409.6) })
        {
            var source = new FakeProcSource().With("/proc/self/status", $"Name:\tx\nVmHWM:\t{anonKb + 1000} kB\nRssAnon:\t{anonKb} kB\nRssFile:\t5 kB\nThreads:\t4\n");
            using var dir = new TempDir();
            LcpChildRun.Outcome r = LcpChildRun.Run(dir, c, LcpChildRun.Options(dir.Sub("scratch"), LcpChildRun.TinyTuning(1), 1, ChildCodec.Packed, c.Sigma) with { ProcSource = source, RssAnonCapMb = cap });
            Assert.That(r.Result.RssVerdict, Is.EqualTo(expected), $"{anonKb} kB against {cap} MB");
            Assert.That(r.Result.RssAnonPeakMb, Is.EqualTo(mb).Within(1e-4), $"{anonKb} kB");
        }
    }

    [Test]
    public void WithRssSamplingOff_TheVerdictIsUnavailable()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(300, 3));
        var source = new FakeProcSource().With("/proc/self/status", "RssAnon:\t1000 kB\n");
        LcpChildRun.Outcome r = LcpChildRun.Run(_dir, c, Tiny(1, sigma: c.Sigma) with { ProcSource = source, MeasureRss = false });
        Assert.That((r.Result.RssAnonPeakMb, r.Result.RssVerdict), Is.EqualTo(((double?)null, RssVerdict.Unavailable)));
        Assert.That(source.Calls, Is.Zero, "nothing was sampled");
    }

    [Test]
    public void OnTheRealLinux_TheRssFigureIsThere_AndTheReservationWasMade()
    {
        if (!OperatingSystem.IsLinux() || !File.Exists("/proc/self/status")) Assert.Ignore("needs /proc");
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(600, 3));
        LcpChildRun.Outcome r = LcpChildRun.Run(_dir, c, Tiny(1, sigma: c.Sigma));
        Assert.That(r.Result.RssAnonPeakMb, Is.GreaterThan(1.0));
        Assert.That(r.Result.RssVerdict, Is.AnyOf(RssVerdict.Within, RssVerdict.Exceeded));
        Assert.That(r.Result.Reservation, Is.EqualTo(ReservationOutcome.Preallocated));
    }

    [Test]
    public void OnASimulatedPlatformWithoutFallocateOrProc_TheStageStillRuns_AndSaysWhatItDidNotCheck()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(600, 3));
        foreach ((long free, ReservationOutcome expected) in new[] { (-1L, ReservationOutcome.Unchecked), (1L << 40, ReservationOutcome.FreeSpaceChecked) })
        {
            using var dir = new TempDir();
            var reservation = new FakeReservationSource { CanPreallocate = false, Free = free };
            LcpChildOptions o = LcpChildRun.Options(dir.Sub("scratch"), LcpChildRun.TinyTuning(2), 2, ChildCodec.Packed, c.Sigma) with { Files = new MappedFileOptions { Source = reservation }, ProcSource = new FakeProcSource() };
            LcpChildRun.Outcome r = LcpChildRun.Run(dir, c, o);
            LcpChildRun.Compare(c, r, ChildCodec.Packed, $"free {free}");
            Assert.That(r.Result.Reservation, Is.EqualTo(expected));
            Assert.That(reservation.PreallocateCalls, Is.Zero);
            Assert.That(reservation.FreeCalls, Is.GreaterThanOrEqualTo(3), "the free space was read for every output");
            Assert.That(r.Result.RssVerdict, Is.EqualTo(RssVerdict.Unavailable));
        }
    }

    [Test]
    public void OnASimulatedPlatform_AFullDiskIsStillAnIoException_FromTheFreeSpaceCheck()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(600, 3));
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        var reservation = new FakeReservationSource { CanPreallocate = false, Free = 100 };
        LcpChildOptions o = Tiny(2, sigma: c.Sigma) with { Files = new MappedFileOptions { Source = reservation } };
        LcpChildRun.AssertFailsCleanly(_dir, sa, textPath, c.N, c.Bits, o, typeof(IOException), "simulated full disk");
    }

    [Test]
    public void TheScalarCommonPrefix_GivesTheSameFiles_AsTheVectorOne()
    {
        string text = StringFamilies.RandomDna(6000, 5) + StringFamilies.RandomDna(900, 12) + StringFamilies.RandomDna(77, 6) + StringFamilies.RandomDna(900, 12);
        LcpChildCase c = LcpChildCase.Of(text);
        string[]? reference = null;
        foreach (bool vectors in new[] { true, false })
        {
            using var dir = new TempDir();
            LcpChildOptions o = LcpChildRun.Options(dir.Sub("scratch"), LcpChildRun.TinyTuning(3) with { LcpChunk = 1000 }, 2, ChildCodec.Packed, c.Sigma);
            LcpChildRun.Outcome r = LcpChildRun.Run(dir, c, o, configure: reader => reader.UseVectors = vectors);
            LcpChildRun.Compare(c, r, ChildCodec.Packed, $"vectors {vectors}");
            string[] h = LcpChildRun.Hashes(r);
            reference ??= h;
            Assert.That(h, Is.EqualTo(reference), $"vectors {vectors}");
        }
    }

    [Test]
    public void TheResult_CarriesTheParameters_TheStats_AndThePeaks()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(3000, 3));
        LcpChildOptions o = Tiny(4, 2, sigma: c.Sigma) with { ArenaGuards = null };
        LcpChildRun.Outcome r = LcpChildRun.Run(_dir, c, o);
        LcpChildResult res = r.Result;
        Assert.That(res.Parameters.N, Is.EqualTo(c.N));
        Assert.That(res.Parameters.LcpChunk, Is.EqualTo(LcpChildRun.TinyTuning(4).LcpChunk));
        Assert.That(res.Lcp.N, Is.EqualTo(c.N));
        Assert.That(res.Scan.N, Is.EqualTo(c.N));
        Assert.That(res.Seconds, Is.GreaterThan(0));
        Assert.That(res.ArenaPeakBytes, Is.GreaterThan(0));
        Assert.That(res.ChildArenaPeakBytes, Is.EqualTo(2L * LcpChildRun.TinyTuning(4).BlockEntries * 16), "the two stack blocks");
        Assert.That(res.Scan.StackPeak, Is.GreaterThan(0));
    }

    [Test]
    public void TheArenaCaps_AreEnforced_PerStage()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(800, 3));
        (string sa, string textPath) = LcpChildRun.WriteInputs(_dir, c);
        LcpChildOptions tooSmallPh4 = Tiny(1, sigma: c.Sigma) with { ArenaCapBytes = 64 };
        LcpChildRun.AssertFailsCleanly(_dir, sa, textPath, c.N, c.Bits, tooSmallPh4, typeof(ArenaBudgetExceededException), "PH4 cap");
        LcpChildOptions tooSmallPh5 = Tiny(1, sigma: c.Sigma) with { ChildArenaCapBytes = 64 };
        LcpChildRun.AssertFailsCleanly(_dir, sa, textPath, c.N, c.Bits, tooSmallPh5, typeof(ArenaBudgetExceededException), "PH5 cap");
    }
}
