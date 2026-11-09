using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Tests.E3Tests;

[TestFixture]
public class FastModTests
{
    [Test]
    public void DivAndMod_EqualTheHardwareOperations_ForEveryDivisorAndEdgeValues()
    {
        var rng = new SplitMix64(5);
        var edges = new uint[] { 0, 1, 2, 3, 1023, 1024, 65535, 65536, int.MaxValue, (uint)int.MaxValue + 1, uint.MaxValue - 1, uint.MaxValue };
        for (uint d = 2; d <= 3000; d++)
        {
            var f = new FastMod32(d);
            foreach (uint a in edges)
            {
                Assert.That(f.Div(a), Is.EqualTo(a / d), $"{a} / {d}");
                Assert.That(f.Mod(a), Is.EqualTo(a % d), $"{a} % {d}");
            }
            for (int i = 0; i < 40; i++)
            {
                uint a = (uint)rng.Next();
                Assert.That(f.Div(a), Is.EqualTo(a / d), $"{a} / {d}");
                Assert.That(f.Mod(a), Is.EqualTo(a % d), $"{a} % {d}");
            }
        }
        foreach (uint d in new uint[] { 4096, 65537, 1u << 20, uint.MaxValue })
        {
            var f = new FastMod32(d);
            for (int i = 0; i < 2000; i++)
            {
                uint a = (uint)rng.Next();
                Assert.That(f.Mod(a), Is.EqualTo(a % d));
                Assert.That(f.Div(a), Is.EqualTo(a / d));
            }
        }
    }

    [Test]
    public void DivisorBelowTwo_IsRejected()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new FastMod32(1));
    }
}

[TestFixture]
public class MappedSinkTests
{
    private string _dir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "emproto-e3-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    [Test]
    public void MappedSaSink_WritesDisjointSlicesFromSeveralThreads_AndTheFileHashEqualsTheMemoryHash()
    {
        const int n = 200_000;
        var data = new uint[n + 1];
        var rng = new SplitMix64(1);
        for (int i = 0; i <= n; i++) data[i] = (uint)rng.Next();
        string path = Path.Combine(_dir, "sa.bin");
        using (var sink = new MappedSaSink(path, n + 1))
        {
            Parallel.For(0, 20, new ParallelOptions { MaxDegreeOfParallelism = 4 }, part =>
            {
                int from = part * 10_000;
                int count = part == 19 ? n + 1 - from : 10_000;
                sink.Write(from, data.AsSpan(from, count));
            });
            Assert.That(sink.MaxWindowBytes, Is.LessThanOrEqualTo(10_001 * 4L));
            Assert.Throws<ArgumentOutOfRangeException>(() => sink.Write(n, new uint[2]));
        }
        byte[] bytes = File.ReadAllBytes(path);
        Assert.That(bytes, Has.Length.EqualTo((n + 1) * 4));
        Assert.That(System.Runtime.InteropServices.MemoryMarshal.Cast<byte, uint>(bytes).ToArray(), Is.EqualTo(data));
        unsafe
        {
            fixed (uint* p = data) Assert.That(SaHash.ComputeFile(path), Is.EqualTo(SaHash.Compute(p + 1, n)));
        }
    }

    [Test]
    public void PositionFiles_AreCountedAndTrackedAsScratch_AndGoAwayWithTheirConsumer()
    {
        ScratchTracker.Reset();
        IoCounters.Reset();
        var w = new PosWriter(_dir, 64, 16);
        for (uint i = 0; i < 10_000; i++) w.Append(i * 3);
        IPosSource src = w.Complete();
        Assert.That(src.Count, Is.EqualTo(10_000));
        Assert.That(Directory.GetFiles(_dir), Has.Length.EqualTo(1));
        Assert.That(ScratchTracker.Live, Is.EqualTo(40_000));
        var buf = new uint[100];
        src.Read(5_000, buf);
        Assert.That(buf[0], Is.EqualTo(15_000u));
        Assert.That(src.At(9_999), Is.EqualTo(29_997u));
        src.Dispose();
        Assert.That(Directory.GetFiles(_dir), Is.Empty);
        Assert.That(ScratchTracker.Live, Is.EqualTo(0));
        Assert.That(ScratchTracker.Peak, Is.EqualTo(40_000));
        Assert.That(IoCounters.WriteBytes, Is.EqualTo(40_000));
    }

    [Test]
    public void AbortedWriter_RemovesItsFile_AndItsTrackedBytes()
    {
        ScratchTracker.Reset();
        var w = new PosWriter(_dir, 64, 16);
        for (uint i = 0; i < 1_000; i++) w.Append(i);
        w.Dispose();
        Assert.That(Directory.GetFiles(_dir), Is.Empty);
        Assert.That(ScratchTracker.Live, Is.EqualTo(0));
    }

    [TestCase("dna:300000", 3000)]
    [TestCase("genome:400000", 2000)]
    [TestCase("tandem:171x400", 500)]
    public void E3Pipeline_PackedFile_RankFile_MappedWindows_MatchesTheInRamKernel(string family, int kDiv)
    {
        if (!File.Exists("/proc/self/stat")) Assert.Ignore("the E3 report reads /proc (Linux)");
        string pk = Path.Combine(_dir, "t.pk");
        string ranks = Path.Combine(_dir, "t.ranks");
        string sa = Path.Combine(_dir, "sa.bin");
        string scratch = Path.Combine(_dir, "scratch");
        var report = new Report(null);
        Assert.That(E3.Pack(Options.Parse(["--family", family, "--out", pk]), report), Is.EqualTo(0));
        Assert.That(E3.MakeRanks(Options.Parse(["--text", pk, "--out", ranks, "--period", "512"]), report), Is.EqualTo(0));
        int rc = E3.Build(Options.Parse(["--text", pk, "--ranks", ranks, "--sa", sa, "--scratch", scratch, "--threads", "3", "--k-div", kDiv.ToString(System.Globalization.CultureInfo.InvariantCulture), "--budget-mb", "16", "--rss-cap-mb", "100000"]), report);
        Assert.That(rc, Is.EqualTo(0), "no scratch file left behind (the RssAnon cap does not apply inside the test host)");

        SuffixTree.EmProto.Gen.SymbolText text = FamilyLoader.Load(family, string.Empty);
        int[] expected = Oracle.SaisKernelSa(text, out _);
        string expectedHash;
        unsafe
        {
            fixed (int* p = expected) expectedHash = SaHash.Compute((uint*)p + 1, expected.Length - 1L);
        }
        Assert.That(SaHash.ComputeFile(sa), Is.EqualTo(expectedHash));
        Assert.That(Directory.GetFiles(scratch), Is.Empty);
    }
}
