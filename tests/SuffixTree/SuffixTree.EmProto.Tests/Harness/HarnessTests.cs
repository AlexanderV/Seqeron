using System.Runtime.InteropServices;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Tests.Harness;

[TestFixture]
public class RssSamplerTests
{
    [Test]
    public unsafe void TouchedNativeMemory_AppearsInRssAnon_AndFreedMemoryDoesNotRaiseTheMaximum()
    {
        if (!RssSampler.IsSupported) Assert.Ignore("needs /proc/self/status with RssAnon (Linux)");
        using var sampler = new RssSampler(5, keepSeries: true);
        sampler.Start();
        RssStats before = sampler.Stop();

        using var second = new RssSampler(5, keepSeries: true);
        second.Start();
        const int mb = 64;
        byte* p = (byte*)NativeMemory.Alloc((nuint)mb << 20);
        for (long i = 0; i < (long)mb << 20; i += 4096) p[i] = 1;
        SpinWait.SpinUntil(static () => false, 60);
        NativeMemory.Free(p);
        RssStats after = second.Stop();

        Assert.That(after.Samples, Is.GreaterThanOrEqualTo(3));
        Assert.That(after.MaxRssAnonKb - before.LastRssAnonKb, Is.GreaterThan(mb * 1024L * 9 / 10));
        Assert.That(second.Series.Count, Is.GreaterThanOrEqualTo(3));
        Assert.That(after.MaxThreads, Is.GreaterThanOrEqualTo(2), "the sampler thread plus the main thread");
    }

    [Test]
    public void StopWithoutStart_AndDoubleStop_AreHarmless()
    {
        var s = new RssSampler(10);
        Assert.That(s.Stop().Samples, Is.EqualTo(0));
        s.Dispose();
    }

    [Test]
    public void StartTwice_Throws()
    {
        if (!RssSampler.IsSupported) Assert.Ignore("needs /proc/self/status with RssAnon (Linux)");
        using var s = new RssSampler(10);
        s.Start();
        Assert.Throws<InvalidOperationException>(s.Start);
    }
}

[TestFixture]
public class IoCountersTests
{
    [Test]
    public void CountedFile_CountsEveryByteOfReadsAndWrites()
    {
        string path = Path.Combine(Path.GetTempPath(), "emproto-io-" + Guid.NewGuid().ToString("N") + ".bin");
        try
        {
            IoCounters.Reset();
            using (var f = new CountedFile(path, FileMode.Create, FileAccess.ReadWrite))
            {
                f.WriteAt(new byte[1000], 0);
                f.WriteAt(new byte[500], 1000);
                Assert.That(f.Length, Is.EqualTo(1500));
                var buf = new byte[700];
                f.ReadExactAt(buf, 0);
                Assert.That(f.ReadAt(buf, 1000), Is.EqualTo(500), "a short read at the end of the file returns the available bytes");
                Assert.Throws<EndOfStreamException>(() => f.ReadExactAt(buf, 1000));
            }
            Assert.That(IoCounters.WriteBytes, Is.EqualTo(1500));
            Assert.That(IoCounters.WriteCalls, Is.EqualTo(2));
            Assert.That(IoCounters.ReadBytes, Is.EqualTo(700 + 500 + 500));
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Test]
    public void ProcIo_ReadsMonotoneCounters()
    {
        if (!File.Exists("/proc/self/io")) Assert.Ignore("needs /proc/self/io (Linux)");
        ProcIoSnapshot a = ProcIoSnapshot.Read();
        string path = Path.Combine(Path.GetTempPath(), "emproto-io-" + Guid.NewGuid().ToString("N") + ".bin");
        try
        {
            File.WriteAllBytes(path, new byte[1 << 20]);
            _ = File.ReadAllBytes(path);
        }
        finally
        {
            File.Delete(path);
        }
        ProcIoSnapshot d = ProcIoSnapshot.Read().Since(a);
        Assert.That(d.WChar, Is.GreaterThanOrEqualTo(1 << 20));
        Assert.That(d.RChar, Is.GreaterThanOrEqualTo(1 << 20));
    }

    [Test]
    public void CountersAgreeWithProcSelfIo_WithinOnePercent()
    {
        if (!File.Exists("/proc/self/io")) Assert.Ignore("needs /proc/self/io (Linux)");
        string path = Path.Combine(Path.GetTempPath(), "emproto-io-" + Guid.NewGuid().ToString("N") + ".bin");
        try
        {
            IoCounters.Reset();
            ProcIoSnapshot before = ProcIoSnapshot.Read();
            using (var f = new CountedFile(path, FileMode.Create, FileAccess.ReadWrite))
            {
                var chunk = new byte[256 * 1024];
                for (int i = 0; i < 48; i++) f.WriteAt(chunk, (long)i * chunk.Length);
                for (int i = 0; i < 48; i++) f.ReadExactAt(chunk, (long)i * chunk.Length);
            }
            ProcIoSnapshot d = ProcIoSnapshot.Read().Since(before);
            long total = 48L * 256 * 1024;
            Assert.That(IoCounters.WriteBytes, Is.EqualTo(total));
            Assert.That(Math.Abs(IoCounters.WriteBytes - d.WChar) / (double)d.WChar, Is.LessThan(0.01));
            Assert.That(Math.Abs(IoCounters.ReadBytes - d.RChar) / (double)d.RChar, Is.LessThan(0.01));
        }
        finally
        {
            File.Delete(path);
        }
    }
}

[TestFixture]
public class ArenaBudgetTests
{
    [Test]
    public void Allocation_IsChargedAndReleased_AndTheCapIsEnforced()
    {
        var budget = new ArenaBudget(1000);
        using var arena = new BudgetedArena(budget);
        _ = arena.Allocate<byte>(400, zero: true, "a");
        int mark = arena.Mark();
        _ = arena.Allocate<byte>(500, zero: true, "b");
        Assert.That(budget.LiveBytes, Is.EqualTo(900));
        Assert.Throws<ArenaBudgetExceededException>(() => arena.Allocate<byte>(101, zero: true, "c"));
        Assert.That(budget.LiveBytes, Is.EqualTo(900), "a refused allocation reserves nothing");
        arena.ReleaseTo(mark);
        Assert.That(budget.LiveBytes, Is.EqualTo(400));
        _ = arena.Allocate<byte>(600, zero: true, "d");
        Assert.That(budget.PeakBytes, Is.EqualTo(1000));
        arena.Dispose();
        Assert.That(budget.LiveBytes, Is.EqualTo(0));
    }

    [Test]
    public void SeveralArenasShareTheBudget_AcrossThreads()
    {
        var budget = new ArenaBudget(10_000);
        int refused = 0;
        Parallel.For(0, 8, worker =>
        {
            _ = worker;
            using var arena = new BudgetedArena(budget);
            for (int i = 0; i < 50; i++)
            {
                try
                {
                    _ = arena.Allocate<byte>(400, zero: false, "x");
                }
                catch (ArenaBudgetExceededException)
                {
                    Interlocked.Increment(ref refused);
                }
            }
        });
        Assert.That(budget.LiveBytes, Is.EqualTo(0));
        Assert.That(budget.PeakBytes, Is.LessThanOrEqualTo(10_000));
        Assert.That(refused, Is.GreaterThan(0), "8 x 50 x 400 bytes is far above the cap");
    }

    [Test]
    public void ZeroLengthBlocks_AreFree()
    {
        var budget = new ArenaBudget(0);
        using var arena = new BudgetedArena(budget);
        _ = arena.Allocate<long>(0, zero: true, "empty");
        Assert.That(budget.LiveBytes, Is.EqualTo(0));
    }
}
