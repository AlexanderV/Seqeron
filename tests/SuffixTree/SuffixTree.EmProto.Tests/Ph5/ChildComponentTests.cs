using System.Runtime.InteropServices;
using System.Text;
using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph5;
using SuffixTree.EmProto.Tests.Build;
using SuffixTree.EmProto.Tests.Ph4;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.Ph5;

/// <summary>
/// Component tests of PH5 written AFTER the code (E5, ADR-004 section 5; the first tests of the gate G1 were written before it, against the stub): the configurations of the scan (stack only, PSV only, packed codec, both), the
/// limits (tagged codec, 32-bit positions, arena budget), the parallel NSV pass, the NEXTSYM pass, the edge cases of the file reader, the stack with another entry type, the backward stack that spills, and the pieces of the
/// harness (the smaps parser, the sampler).
/// </summary>
[TestFixture]
public class ChildComponentTests
{
    private string _dir = "";

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "e5c-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    private string Scratch => Path.Combine(_dir, "scratch");

    [Test]
    public void Configurations_StackOnly_PsvOnly_PackedPlusPsv_GiveTheSameStatisticsAndTheirOwnFiles()
    {
        ChildReference r = ChildReference.Of(StringFamilies.Fibonacci(4000));
        (string b, string o) = r.WriteLcpFiles(_dir, "lcp");
        string child = Path.Combine(_dir, "c.bin"), psv = Path.Combine(_dir, "p.bin");
        ChildScanResult stackOnly = ExternalChild.BuildFromFiles(r.N, b, o, new ChildBuildOutputs(), ChildChecks.Tiny(Scratch)).Forward;
        ChildChecks.CompareStats(r, stackOnly, "stack only");
        Assert.That(stackOnly.ChildWindows, Is.EqualTo(default(WindowStats)), "no writer without an output");
        ChildScanResult psvOnly = ExternalChild.BuildFromFiles(r.N, b, o, new ChildBuildOutputs { Psv = psv }, ChildChecks.Tiny(Scratch)).Forward;
        ChildChecks.CompareStats(r, psvOnly, "psv only");
        ChildChecks.CompareArrays(r.Psv, ChildFiles.ReadSmaller(psv), "PSV", "psv only");
        Assert.That(psvOnly.ChildWindows, Is.EqualTo(default(WindowStats)));
        Assert.That(psvOnly.PsvWindows.Writes, Is.EqualTo(r.N + 1));
        ChildScanResult packedPsv = ExternalChild.BuildFromFiles(r.N, b, o, new ChildBuildOutputs { Child = child, Psv = psv }, ChildChecks.Tiny(Scratch, codec: ChildCodec.Packed)).Forward;
        ChildChecks.CompareChild(r, child, ChildCodec.Packed, "packed + psv");
        ChildChecks.CompareArrays(r.Psv, ChildFiles.ReadSmaller(psv), "PSV", "packed + psv");
        Assert.That(packedPsv.ChildWindows.Writes, Is.LessThan(r.N + 1), "the packed codec writes fewer fields (a next link only when it is not k + 1)");
        Assert.That(packedPsv.PsvWindows.Writes, Is.EqualTo(r.N + 1), "one PSV entry per rank");
    }

    [Test]
    public void ParallelBackwardPass_GivesTheSameFiles_AndReportsBothPasses()
    {
        foreach (string text in new[] { StringFamilies.AllSame(900) + "b", StringFamilies.Fibonacci(3000), StringFamilies.RandomDna(2000, 1), ChildSamples.U8WithEscapes() })
        {
            ChildReference r = ChildReference.Of(text);
            ChildBuildStats seq = ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch, blockEntries: 5), _dir, "sequential");
            string nsv = File.ReadAllText(Path.Combine(_dir, "nsv.bin"), Encoding.Latin1);
            (string b, string o) = r.WriteLcpFiles(_dir, "lcp");
            string nsv2 = Path.Combine(_dir, "nsv2.bin");
            ChildBuildStats par = ExternalChild.BuildFromFiles(r.N, b, o, new ChildBuildOutputs { Nsv = nsv2 }, ChildChecks.Tiny(Scratch, blockEntries: 5) with { ParallelBackward = true });
            Assert.That(File.ReadAllText(nsv2, Encoding.Latin1), Is.EqualTo(nsv), "the concurrent pass writes the same NSV");
            Assert.That(par.Backward, Is.Not.Null);
            Assert.That(par.BackwardConcurrent, Is.True);
            Assert.That(seq.BackwardConcurrent, Is.False);
            Assert.That(par.LcpBytes, Is.EqualTo(new FileInfo(b).Length + new FileInfo(o).Length), "the volume each pass reads");
            Assert.That(par.Backward!.StackPeak, Is.EqualTo(seq.Backward!.StackPeak));
            Assert.That(par.Forward.NodeCount, Is.EqualTo(seq.Forward.NodeCount));
        }
    }

    [Test]
    public void ParallelBackwardPass_DoesNotHang_WhenTheForwardPassFails()
    {
        ChildReference r = ChildReference.Of(ChildSamples.U8WithEscapes());
        (string b, string o) = r.WriteLcpFiles(_dir, "lcp");
        byte[] pairs = File.ReadAllBytes(o);
        pairs[0] ^= 1;
        string bad = Path.Combine(_dir, "bad.ovf");
        File.WriteAllBytes(bad, pairs);
        Assert.Throws<InvalidDataException>(() => ExternalChild.BuildFromFiles(r.N, b, bad, new ChildBuildOutputs { Child = Path.Combine(_dir, "c"), Nsv = Path.Combine(_dir, "n") }, ChildChecks.Tiny(Scratch) with { ParallelBackward = true }));
    }

    [Test]
    public void BackwardStack_SpillsWhenTheLcpDecreasesToTheLeft_AndTheNsvStillEqualsTheOracle()
    {
        // a^k b: the suffixes in order are the empty one, a^k b, a^(k-1) b, ..., b: the LCP falls from k - 1 to 0, so the backward scan (from the right) climbs 0, 1, ..., k - 1
        ChildReference r = ChildReference.Of(new string('a', 700) + "b");
        foreach (int block in new[] { 4, 7, 64 })
        {
            ChildBuildStats st = ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch, blockEntries: block), _dir, $"a^k b, block {block}");
            NsvResult back = st.Backward!;
            Assert.That(back.StackPeak, Is.EqualTo(700), "height = max LCP + 1");
            long blocks = Math.Max(0, ((700 - 1) / block) - 1);
            Assert.That(back.Spill.BlocksWritten, Is.EqualTo(blocks));
            Assert.That(back.Spill.BytesWritten, Is.EqualTo(blocks * block * 8), "8-byte entries");
            Assert.That(back.Spill.BlocksRead, Is.EqualTo(blocks));
            Assert.That(st.Forward.StackPeak, Is.LessThan(5), "the forward stack stays small on this text");
            Assert.That(back.Windows.Flushes, Is.GreaterThan(0), "the NSV windows are flushed when released");
            Assert.That(back.Windows.Writes, Is.EqualTo(r.N + 1), "one NSV entry per rank");
        }
    }

    [Test]
    public void TheSpillStack_WorksForAnyUnmanagedEntry_AndEveryBlockIsWrittenAtItsOwnOffset()
    {
        using var arena = new BudgetedArena(new ArenaBudget(long.MaxValue), ArenaOptions.Guarded(0xCD));
        string path = Path.Combine(_dir, "v.spill");
        using var stack = new SpillStack<ValueEntry>(arena, path, 5);
        for (uint i = 0; i < 53; i++) stack.Push(new ValueEntry { Value = i, Index = i * 7 });
        Assert.That(new FileInfo(path).Length, Is.EqualTo(stack.Stats.PeakFileBytes));
        Assert.That(stack.Stats.PeakFileBytes % (5 * 8), Is.Zero, "whole blocks of 5 entries of 8 bytes");
        byte[] file = File.ReadAllBytes(path);
        ReadOnlySpan<ValueEntry> onDisk = MemoryMarshal.Cast<byte, ValueEntry>(file);
        for (int i = 0; i < onDisk.Length; i++) Assert.That((onDisk[i].Value, onDisk[i].Index), Is.EqualTo(((uint)i, (uint)i * 7)), $"entry {i} of the file");
    }

    [Test]
    public void ScanLimits_TaggedCodecAndThe32BitPositions_AreRefusedBeforeAnythingIsCreated()
    {
        string c = Path.Combine(_dir, "c.bin");
        Assert.Throws<NotSupportedException>(() => ChildScan.Create(1L << 30, new ChildScanOutputs { ChildPath = c }, ChildChecks.Tiny(Scratch)));
        Assert.That(File.Exists(c), Is.False);
        Assert.Throws<NotSupportedException>(() => ChildScan.Create(uint.MaxValue - 1L, new ChildScanOutputs(), ChildChecks.Tiny(Scratch, codec: ChildCodec.Packed)));
        Assert.Throws<ArgumentOutOfRangeException>(() => ChildScan.Create(-1, new ChildScanOutputs(), ChildChecks.Tiny(Scratch)));
        Assert.Throws<ArgumentNullException>(() => ChildScan.Create(5, null!, ChildChecks.Tiny(Scratch)));
    }

    [Test]
    public void TheArenaBudget_CapsTheTwoBlocks_AndEverythingIsReleasedWhenTheCapFires()
    {
        var budget = new ArenaBudget(1000);
        ChildBuildOptions o = ChildChecks.Tiny(Scratch, blockEntries: 100) with { Budget = budget };
        Assert.Throws<ArenaBudgetExceededException>(() => ChildScan.Create(10, new ChildScanOutputs { ChildPath = Path.Combine(_dir, "c.bin") }, o));
        Assert.That(budget.LiveBytes, Is.Zero, "the arena is disposed, the budget balanced");
        var enough = new ArenaBudget(2 * 100 * 16);
        using (ChildScan.Create(10, new ChildScanOutputs(), ChildChecks.Tiny(Scratch, blockEntries: 100) with { Budget = enough }))
        {
            Assert.That(enough.LiveBytes, Is.EqualTo(3200), "two blocks of 100 entries of 16 bytes, nothing else");
        }
        Assert.That(enough.LiveBytes, Is.Zero);
    }

    [Test]
    public void EmptyAndTinyTexts_AreHandled()
    {
        foreach (string text in new[] { "a", "ab", "aa" })
        {
            ChildReference r = ChildReference.Of(text);
            ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch), _dir, text);
        }
        // the empty text: n = 0 has one rank (the empty suffix), no interval, NODE_COUNT = 1
        string child = Path.Combine(_dir, "e.child");
        using ChildScan scan = ChildScan.Create(0, new ChildScanOutputs { ChildPath = child }, ChildChecks.Tiny(Scratch));
        scan.Consume(0, new uint[] { 0 });
        ChildScanResult res = scan.Complete();
        Assert.That((res.NodeCount, res.InternalCount, res.StackPeak), Is.EqualTo((1L, 0L, 1L)));
        Assert.That(File.ReadAllBytes(child), Is.EqualTo(new byte[4]), "one zero entry, as the in-RAM layout of the empty text");
    }

    [Test]
    public void CompleteTwice_IsRefused_AndDisposeAfterCompleteIsHarmless()
    {
        ChildReference r = ChildReference.Of("banana");
        uint[] dense = ChildFiles.Dense(r.Lcp);
        ChildScan scan = ChildScan.Create(r.N, new ChildScanOutputs(), ChildChecks.Tiny(Scratch));
        scan.Consume(0, dense);
        _ = scan.Complete();
        Assert.Throws<InvalidOperationException>(() => scan.Complete());
        scan.Dispose();
        scan.Dispose();
        Assert.Throws<ObjectDisposedException>(() => scan.Consume(0, dense));
    }

    [Test]
    public void ValuesAboveN_AndABadFirstValue_AreRefusedByTheScan()
    {
        using ChildScan scan = ChildScan.Create(3, new ChildScanOutputs(), ChildChecks.Tiny(Scratch));
        Assert.Throws<InvalidDataException>(() => scan.Consume(0, new uint[] { 1 }), "LCP[0] must be 0");
        using ChildScan scan2 = ChildScan.Create(3, new ChildScanOutputs(), ChildChecks.Tiny(Scratch));
        Assert.Throws<InvalidDataException>(() => scan2.Consume(0, new uint[] { 0, 0, 9 }), "a value above n");
        using ChildScan scan3 = ChildScan.Create(3, new ChildScanOutputs(), ChildChecks.Tiny(Scratch));
        Assert.Throws<InvalidDataException>(() => scan3.Consume(0, new uint[] { 0, 0, 4 }), "n + 1 is above n");
        using ChildScan scan4 = ChildScan.Create(3, new ChildScanOutputs(), ChildChecks.Tiny(Scratch));
        Assert.DoesNotThrow(() => scan4.Consume(0, new uint[] { 0, 0, 3 }), "n itself is not refused by the scan");
    }

    [Test]
    public void TaggedCodec_AcceptsTheLargestTextItCanHold()
    {
        // 2^30 - 1 symbols: the file is 4 GiB sparse (never touched); the limit is exact
        string c = Path.Combine(_dir, "big.cld");
        using ChildScan scan = ChildScan.Create((1L << 30) - 1, new ChildScanOutputs { ChildPath = c }, ChildChecks.Tiny(Scratch));
        Assert.That(new FileInfo(c).Length, Is.EqualTo(4L * (1L << 30)));
    }

    [Test]
    public void CompactLcpFile_RefusesWhatDoesNotFitItsBaseArray()
    {
        string b = Path.Combine(_dir, "b.lcp"), o = Path.Combine(_dir, "b.ovf");
        File.WriteAllBytes(b, new byte[4 * 6]);
        File.WriteAllBytes(o, new byte[8]);
        Assert.Throws<ArgumentException>(() => new CompactLcpFile(b, o, 5), "a 4-byte base array has no overflow list");
        foreach (int length in new[] { 5, 7, 11, 13, 25 })
        {
            File.WriteAllBytes(b, new byte[length]);
            File.WriteAllBytes(o, []);
            Assert.Throws<ArgumentException>(() => new CompactLcpFile(b, o, 5), $"a base file of {length} bytes is not 6 entries of 1, 2 or 4 bytes");
        }
        File.WriteAllBytes(b, new byte[4 * 6]);
        File.WriteAllBytes(o, new byte[8]);
        File.WriteAllBytes(b, new byte[6]);
        foreach (int bad in new[] { 4, 5, 12 })
        {
            File.WriteAllBytes(o, new byte[bad]);
            Assert.Throws<ArgumentException>(() => new CompactLcpFile(b, o, 5), $"an overflow file of {bad} bytes is not a whole number of pairs");
        }
        File.WriteAllBytes(o, []);
        using var ok = new CompactLcpFile(b, o, 5);
        Assert.That(ok.Width, Is.EqualTo(1));
        // a pair below the escape value, and a pair that no escape slot uses
        File.WriteAllBytes(b, [0, 0, 255, 0, 0, 0]);
        File.WriteAllBytes(o, MemoryMarshal.AsBytes(new uint[] { 2, 100 }.AsSpan()).ToArray());
        using var low = new CompactLcpFile(b, o, 5);
        Assert.Throws<InvalidDataException>(() => low.ForEachForward((_, _) => { }), "the value of the pair is below 255");
        Assert.Throws<InvalidDataException>(() => low.ForEachBackward((_, _) => { }));
        File.WriteAllBytes(b, [0, 0, 0, 0, 0, 0]);
        File.WriteAllBytes(o, MemoryMarshal.AsBytes(new uint[] { 2, 300 }.AsSpan()).ToArray());
        using var extra = new CompactLcpFile(b, o, 5);
        Assert.Throws<InvalidDataException>(() => extra.ForEachForward((_, _) => { }), "a pair without an escape");
        Assert.Throws<InvalidDataException>(() => extra.ForEachBackward((_, _) => { }));
    }

    [Test]
    public void CompactLcpFile_ReadsTheOverflowInWindows_ForwardAndBackward()
    {
        // 700 escapes read through a 16-byte buffer (2 pairs per window): the windows move in both directions
        ChildReference r = ChildReference.Of(ChildSamples.U8WithEscapes());
        (string b, string o) = r.WriteLcpFiles(_dir, "w");
        Assert.That(r.Lcp.OverflowCount, Is.GreaterThan(100));
        using var file = new CompactLcpFile(b, o, r.N, 16);
        uint[] dense = ChildFiles.Dense(r.Lcp);
        var seen = new List<uint>();
        long calls0 = IoCounters.ReadCalls;
        file.ForEachForward((_, v) => seen.AddRange(v.ToArray()));
        long forward = IoCounters.ReadCalls - calls0;
        Assert.That(seen, Is.EqualTo(dense));
        seen.Clear();
        calls0 = IoCounters.ReadCalls;
        file.ForEachBackward((_, v) => seen.InsertRange(0, v.ToArray()));
        long backward = IoCounters.ReadCalls - calls0;
        Assert.That(seen, Is.EqualTo(dense));
        // the overflow is read in windows of two pairs that move with the cursor: about one read per two escapes in either direction (a window that did not follow the direction would read once per escape)
        Assert.That(backward, Is.LessThanOrEqualTo(forward + 4), $"backward {backward} reads, forward {forward}");
    }

    [Test]
    public void NextSym_IsTheSymbolAfterTheCommonPrefix_AndTheEndMarkAtTheEndOfTheText()
    {
        foreach (string text in new[] { "banana", "mississippi", StringFamilies.Fibonacci(2000), StringFamilies.RandomDna(3000, 2), new string('a', 300), StringFamilies.Runs(2500, 16, 4) })
        {
            ChildReference r = ChildReference.Of(text);
            DenseCase c = BuilderTestSupport.Dense("t", text);
            using PackedText packed = PackedText.Pack(c.Symbols, c.MinBits);
            (string b, string o) = r.WriteLcpFiles(_dir, "ns");
            string sa = LcpHarness.WriteSa(r.Index.SuffixArray.Select(x => (uint)x).ToArray(), _dir);
            string outPath = Path.Combine(_dir, "ns.out");
            NextSymResult res;
            using (var file = new CompactLcpFile(b, o, r.N, 64)) res = NextSymScan.Run(packed, sa, file, outPath, ChildChecks.Tiny(Scratch));
            byte[] got = File.ReadAllBytes(outPath);
            int[] lcp = r.Lcp.ToArray();
            Assert.That(got.Length, Is.EqualTo(r.N + 1));
            for (int i = 0; i <= r.N; i++)
            {
                int p = r.Index.SuffixArray[i] + lcp[i];
                byte expected = p < r.N ? c.Symbols[p] : (byte)0xFF;
                Assert.That(got[i], Is.EqualTo(expected), $"rank {i}");
            }
            Assert.That(got[0], Is.EqualTo(NextSymScan.NoSymbol), "rank 0 is the empty suffix: its position is n");
            Assert.That(res.TextReads, Is.EqualTo(r.N + 1));
        }
    }

    [Test]
    public void NextSym_RefusesAWrongSuffixArrayFileAndAWideAlphabet()
    {
        ChildReference r = ChildReference.Of("banana");
        DenseCase c = BuilderTestSupport.Dense("t", "banana");
        using PackedText packed = PackedText.Pack(c.Symbols, c.MinBits);
        (string b, string o) = r.WriteLcpFiles(_dir, "ns");
        string sa = Path.Combine(_dir, "short.sa");
        File.WriteAllBytes(sa, new byte[12]);
        using var file = new CompactLcpFile(b, o, r.N);
        Assert.Throws<ArgumentException>(() => NextSymScan.Run(packed, sa, file, Path.Combine(_dir, "x"), ChildChecks.Tiny(Scratch)));
        using PackedText wide = PackedText.Pack(new int[] { 1, 2, 300, 4 }, 9);
        string saWide = Path.Combine(_dir, "wide.sa");
        File.WriteAllBytes(saWide, new byte[5 * 4]); // the right length for n = 4: only the width of the alphabet is wrong
        ArgumentException ex = Assert.Throws<ArgumentException>(() => NextSymScan.Run(wide, saWide, file, Path.Combine(_dir, "y"), ChildChecks.Tiny(Scratch)))!;
        Assert.That(ex.Message, Does.Contain("8 bits"));
    }

    [Test]
    public void SmapsParser_SumsTheDirtyPagesOfTheMappingsUnderThePrefixOnly()
    {
        string text = """
            7f0000000000-7f0000100000 rw-s 00000000 fe:00 123      /tmp/run/out.cld
            Size:               1024 kB
            Private_Clean:         0 kB
            Private_Dirty:       300 kB
            Shared_Dirty:         12 kB
            Rss:                 312 kB
            7f0000200000-7f0000300000 rw-s 00100000 fe:00 123      /tmp/run/out.cld
            Private_Dirty:       100 kB
            Shared_Dirty:          0 kB
            7f0000400000-7f0000500000 rw-s 00000000 fe:00 124      /tmp/run/other.bin
            Private_Dirty:       999 kB
            Shared_Dirty:        999 kB
            7f0000600000-7f0000700000 rw-p 00000000 00:00 0
            Private_Dirty:       555 kB
            7ffd0000000-7ffd0100000 rw-p 00000000 00:00 0          [stack]
            Private_Dirty:       777 kB
            """;
        Assert.That(DirtySampler.SumFileMapped(Encoding.UTF8.GetBytes(text), Encoding.UTF8.GetBytes("/tmp/run/out")), Is.EqualTo(412));
        Assert.That(DirtySampler.SumFileMapped(Encoding.UTF8.GetBytes(text), Encoding.UTF8.GetBytes("/tmp/run/other")), Is.EqualTo(1998));
        Assert.That(DirtySampler.SumFileMapped(Encoding.UTF8.GetBytes(text), Encoding.UTF8.GetBytes("/nowhere")), Is.Zero);
        Assert.That(DirtySampler.SumFileMapped([], Encoding.UTF8.GetBytes("/x")), Is.Zero);
    }

    [Test]
    public void DirtySampler_TakesABaseline_SeesDirtyPagesOfAMappedFile_AndMarks()
    {
        if (!File.Exists("/proc/meminfo")) Assert.Ignore("Linux only");
        string path = Path.Combine(_dir, "mapped.bin");
        using var sampler = new DirtySampler(5, filePrefix: path, keepSeries: true);
        sampler.Start();
        sampler.Mark();
        using (var w = new AlignedWindowWriter(path, 16 << 20, 1 << 20, 16, flushOnRelease: false))
        {
            for (long i = 0; i < (16 << 20) / 4; i += 1024) w.Write32(i, 1);
            SpinWait.SpinUntil(() => false, 80); // let the sampler see the dirty pages
            DirtyStats during = sampler.Stop();
            Assert.That(during.Samples, Is.GreaterThan(3));
            Assert.That(during.MaxProcessDirtyKb, Is.GreaterThan(8 * 1024), "the per-process view sees the mapped dirty pages of the file under the prefix");
            Assert.That(during.MaxDeltaKb, Is.GreaterThanOrEqualTo(0));
            Assert.That(during.MarkBaselineKb, Is.GreaterThanOrEqualTo(0));
        }
        Assert.That(sampler.Series, Is.Not.Empty);
    }

    [Test]
    public void Healthy_PinsEveryTerm()
    {
        Assert.That(E5.Healthy(10, 256, 0, 0, 0), Is.True);
        Assert.That(E5.Healthy(300, 256, 0, 0, 0), Is.False, "RssAnon above the cap");
        Assert.That(E5.Healthy(256, 256, 0, 0, 0), Is.True, "at the cap");
        Assert.That(E5.Healthy(10, 256, 1, 0, 0), Is.False, "a scratch entry left");
        Assert.That(E5.Healthy(10, 256, 0, 1, 0), Is.False, "tracked scratch bytes left");
        Assert.That(E5.Healthy(10, 256, 0, 0, 1), Is.False, "arena bytes left");
    }

    private (string Text, string Sa, string Lcp) Prepare(string family)
    {
        string pk = Path.Combine(_dir, "t.pk");
        string ranks = Path.Combine(_dir, "t.ranks");
        string sa = Path.Combine(_dir, "t.sa");
        var report = new Report(null);
        Assert.That(E3.Pack(Options.Parse(["--family", family, "--out", pk]), report), Is.EqualTo(0));
        Assert.That(E3.MakeRanks(Options.Parse(["--text", pk, "--out", ranks, "--period", "512"]), report), Is.EqualTo(0));
        Assert.That(E3.Build(Options.Parse(["--text", pk, "--ranks", ranks, "--sa", sa, "--scratch", Path.Combine(_dir, "e3scratch"), "--threads", "2", "--k-div", "500", "--budget-mb", "16", "--rss-cap-mb", "100000"]), report), Is.EqualTo(0));
        string lcp = Path.Combine(_dir, "t.oracle");
        Assert.That(E4.Oracle(Options.Parse(["--text", pk, "--sa", sa, "--out", lcp]), report), Is.EqualTo(0));
        return (pk, sa, lcp);
    }

    [TestCase("genome:200000", "a", "tagged", false)]
    [TestCase("genome:200000", "f", "tagged", true)]
    [TestCase("genome:200000", "both", "packed", false)]
    [TestCase("a:70000", "both", "tagged", true)]
    [TestCase("a:70000", "a", "packed", false)]
    [TestCase("fib:150000", "both", "tagged", false)]
    [TestCase("tandem:171x300", "f", "tagged", true)]
    [TestCase("nrun:100000:60000", "both", "packed", true)]
    [TestCase("genome:100000", "stats", "tagged", false)]
    public void Pipeline_Lcp_Oracle_E5_Verify_Nextsym(string family, string layout, string codec, bool parallel)
    {
        if (!File.Exists("/proc/self/stat")) Assert.Ignore("the E5 report reads /proc (Linux)");
        (string pk, string sa, string lcp) = Prepare(family);
        string rows = Path.Combine(_dir, "rows.tsv");
        var report = new Report(rows);
        string reference = Path.Combine(_dir, "ref");
        Assert.That(E5.Oracle(Options.Parse(["--text", pk, "--lcp", lcp, "--out", reference, "--sa", sa]), report), Is.EqualTo(0));
        List<string> args = ["--text", pk, "--lcp", lcp, "--out", Path.Combine(_dir, "out"), "--scratch", Path.Combine(_dir, "scratch"), "--layout", layout, "--codec", codec, "--verify", reference, "--case", "t", "--block-entries", "64", "--window-kb", "4", "--slots", "3", "--io-kb", "2", "--spill-flush-mb", "0", "--sa", sa, "--nextsym", "--rss-cap-mb", "100000"];
        if (parallel) args.Add("--parallel-nsv");
        Assert.That(E5.Build(Options.Parse([.. args]), report), Is.EqualTo(0), "equal to the in-RAM reference, no scratch left, budget balanced");
        string text = File.ReadAllText(rows);
        if (layout is "a" or "both") Assert.That(text, Does.Contain("\tchild_equals_oracle\tyes"));
        if (layout is "f" or "both") Assert.That(text, Does.Contain("\tpsv_equals_oracle\tyes").And.Contain("\tnsv_equals_oracle\tyes"));
        Assert.That(text, Does.Contain("\tstats_equal_oracle\tyes").And.Contain("\tnextsym_equals_oracle\tyes").And.Contain("\tspill_expected\t"));
        Assert.That(text, Does.Contain($"parallel NSV={parallel}"));
        if (layout != "stats") Assert.That(text, Does.Not.Contain(" 0 flushes covering 0 bytes"), "the windows are flushed when released");
        Assert.That(Directory.GetFileSystemEntries(Path.Combine(_dir, "scratch")), Is.Empty);
    }

    [Test]
    public void NoFlushOption_TurnsTheFlushOfReleasedWindowsOff()
    {
        if (!File.Exists("/proc/self/stat")) Assert.Ignore("the E5 report reads /proc (Linux)");
        (string pk, _, string lcp) = Prepare("genome:100000");
        string rows = Path.Combine(_dir, "rows.tsv");
        var report = new Report(rows);
        string[] args = ["--text", pk, "--lcp", lcp, "--out", Path.Combine(_dir, "out"), "--scratch", Path.Combine(_dir, "scratch"), "--layout", "a", "--window-kb", "4", "--slots", "2", "--no-flush", "--rss-cap-mb", "100000", "--case", "noflush"];
        Assert.That(E5.Build(Options.Parse(args), report), Is.EqualTo(0));
        Assert.That(File.ReadAllText(rows), Does.Contain(" 0 flushes covering 0 bytes"));
        string withFlush = Path.Combine(_dir, "rows2.tsv");
        Assert.That(E5.Build(Options.Parse([.. args.Where(a => a != "--no-flush"), "--case", "flush"]), new Report(withFlush)), Is.EqualTo(0));
        Assert.That(File.ReadAllText(withFlush), Does.Not.Contain(" 0 flushes covering 0 bytes"));
    }

    [TestCase("genome:100000", "a")]
    [TestCase("a:60000", "f")]
    [TestCase("tandem:171x200", "both")]
    public void Fused_Ph4Hook_GivesTheSameFilesAsTheReference(string family, string layout)
    {
        if (!File.Exists("/proc/self/stat")) Assert.Ignore("the E5 report reads /proc (Linux)");
        (string pk, string sa, string lcp) = Prepare(family);
        string rows = Path.Combine(_dir, "rows.tsv");
        var report = new Report(rows);
        string reference = Path.Combine(_dir, "ref");
        Assert.That(E5.Oracle(Options.Parse(["--text", pk, "--lcp", lcp, "--out", reference]), report), Is.EqualTo(0));
        string[] args = ["--fused", "--text", pk, "--sa", sa, "--out", Path.Combine(_dir, "out"), "--scratch", Path.Combine(_dir, "scratch"), "--scale-n", "3000000000", "--threads", "3", "--layout", layout, "--verify", reference, "--case", "t", "--block-entries", "32", "--window-kb", "4", "--rss-cap-mb", "100000"];
        Assert.That(E5.Build(Options.Parse(args), report), Is.EqualTo(0));
        string text = File.ReadAllText(rows);
        Assert.That(text, Does.Contain("fused with PH4=True").And.Contain("\tph4_parameters\t").And.Contain("\tstats_equal_oracle\tyes"));
        Assert.That(Directory.GetFileSystemEntries(Path.Combine(_dir, "scratch")), Is.Empty, "the PH4 directory is removed once it is empty");
    }

    [Test]
    public void Verify_AgainstAWrongOracle_ReturnsOne_AndSaysNo()
    {
        if (!File.Exists("/proc/self/stat")) Assert.Ignore("the E5 report reads /proc (Linux)");
        (string pk, _, string lcp) = Prepare("genome:100000");
        string rows = Path.Combine(_dir, "rows.tsv");
        var report = new Report(rows);
        string reference = Path.Combine(_dir, "ref");
        Assert.That(E5.Oracle(Options.Parse(["--text", pk, "--lcp", lcp, "--out", reference]), report), Is.EqualTo(0));
        string[] correctStats = File.ReadAllLines(reference + ".stats");
        byte[] child = File.ReadAllBytes(reference + ".cldT");
        child[4001] ^= 1;
        File.WriteAllBytes(reference + ".cldT", child);
        string[] args = ["--text", pk, "--lcp", lcp, "--out", Path.Combine(_dir, "out"), "--scratch", Path.Combine(_dir, "scratch"), "--layout", "a", "--verify", reference, "--rss-cap-mb", "100000"];
        Assert.That(E5.Build(Options.Parse(args), report), Is.EqualTo(1));
        Assert.That(File.ReadAllText(rows), Does.Contain("\tchild_equals_oracle\tNO\t").And.Contain("first difference at byte 4001"));
        File.WriteAllLines(reference + ".stats", File.ReadAllLines(reference + ".stats").Select(l => l.StartsWith("nodecount=", StringComparison.Ordinal) ? "nodecount=1" : l));
        File.WriteAllBytes(reference + ".cldT", MemoryMarshal.AsBytes(ChildFiles.ReadUInts(Path.Combine(_dir, "out.cld")).AsSpan()).ToArray());
        Assert.That(E5.Build(Options.Parse(args), report), Is.EqualTo(1), "the statistics differ");
        Assert.That(File.ReadAllText(rows), Does.Contain("\tstats_equal_oracle\tNO\t"));
        // each statistic on its own: NODE_COUNT, the internal-node count, LRS_DEPTH and LRS_RANK each decide the comparison
        foreach (string key in new[] { "nodecount", "internal", "lrsdepth", "lrsrank" })
        {
            File.WriteAllLines(reference + ".stats", correctStats.Select(l => l.StartsWith(key + "=", StringComparison.Ordinal) ? key + "=12345678" : l));
            Assert.That(E5.Build(Options.Parse(args), report), Is.EqualTo(1), $"a wrong {key} fails the run");
        }
        File.WriteAllLines(reference + ".stats", correctStats);
        Assert.That(E5.Build(Options.Parse(args), report), Is.EqualTo(0), "with the right statistics and the output as reference the run passes");
        Assert.That(E5.Build(Options.Parse([.. args[..^2], "--rss-cap-mb", "0.001"]), report), Is.EqualTo(1), "the RssAnon cap of the gate");
    }

    [Test]
    public void NextSym_AgainstAWrongOracle_FailsTheRun()
    {
        if (!File.Exists("/proc/self/stat")) Assert.Ignore("the E5 report reads /proc (Linux)");
        (string pk, string sa, string lcp) = Prepare("genome:100000");
        string rows = Path.Combine(_dir, "rows.tsv");
        var report = new Report(rows);
        string reference = Path.Combine(_dir, "ref");
        Assert.That(E5.Oracle(Options.Parse(["--text", pk, "--lcp", lcp, "--out", reference, "--sa", sa]), report), Is.EqualTo(0));
        byte[] nextSym = File.ReadAllBytes(reference + ".nsym");
        nextSym[777] ^= 1;
        File.WriteAllBytes(reference + ".nsym", nextSym);
        string[] args = ["--text", pk, "--lcp", lcp, "--out", Path.Combine(_dir, "out"), "--scratch", Path.Combine(_dir, "scratch"), "--layout", "stats", "--nextsym", "--sa", sa, "--verify", reference, "--rss-cap-mb", "100000"];
        Assert.That(E5.Build(Options.Parse(args), report), Is.EqualTo(1));
        Assert.That(File.ReadAllText(rows), Does.Contain("\tnextsym_equals_oracle\tNO\t").And.Contain("first difference at byte 777"));
    }

    [Test]
    public void Commands_RejectBadOptions_AndTheCliDispatchesThem()
    {
        var report = new Report(null);
        Assert.Throws<ArgumentException>(() => E5.Build(Options.Parse(["--text", "x", "--lcp", "y", "--out", "z", "--scratch", "s", "--layout", "z"]), report));
        Assert.Throws<ArgumentException>(() => E5.Build(Options.Parse(["--text", "x", "--lcp", "y", "--out", "z", "--scratch", "s", "--codec", "wide"]), report));
        Assert.Throws<ArgumentException>(() => E5.Build(Options.Parse(["--text", "x"]), report));
        Assert.That(Cli.Run(["e5", "--text", "x"]), Is.EqualTo(2));
        Assert.That(Cli.Run(["e5-oracle"]), Is.EqualTo(2));
        Assert.That(E5.ChildFile("p"), Is.EqualTo("p.cld"));
        Assert.That(E5.PsvFile("p"), Is.EqualTo("p.psv"));
        Assert.That(E5.NsvFile("p"), Is.EqualTo("p.nsv"));
        var sw = new StringWriter();
        TextWriter old = Console.Out;
        Console.SetOut(sw);
        try
        {
            Assert.That(Cli.Run(["--help"]), Is.EqualTo(0));
        }
        finally
        {
            Console.SetOut(old);
        }
        Assert.That(sw.ToString(), Does.Contain("e5 ").And.Contain("e5-oracle").And.Contain("--parallel-nsv").And.Contain("--nextsym"));
    }
}
