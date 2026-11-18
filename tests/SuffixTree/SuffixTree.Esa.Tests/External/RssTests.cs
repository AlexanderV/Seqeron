using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>RssAnon through an injectable source (G4: a platform without /proc reports "unavailable" and never throws).</summary>
[TestFixture]
public class RssReaderTests
{
    private const string Status = "Name:\tdotnet\nVmPeak:\t 999 kB\nVmHWM:\t  4321 kB\nRssAnon:\t   1234 kB\nRssFile:\t    567 kB\nRssShmem:\t 0 kB\nThreads:\t12\n";

    [Test]
    public void ReadsTheFields_OfARealisticStatusFile()
    {
        var r = new RssReader(new FakeProcSource().With("/proc/self/status", Status));
        RssSnapshot s = r.Read();
        Assert.That(s, Is.EqualTo(new RssSnapshot(true, 1234, 567, 4321, 12)));
        Assert.That(s.RssAnonMb, Is.EqualTo(1234 * 1024.0 / 1e6).Within(1e-12));
    }

    [Test]
    public void NoProc_IsUnavailable()
    {
        var r = new RssReader(new FakeProcSource());
        RssSnapshot s = r.Read();
        Assert.That(s.Available, Is.False);
        Assert.That(s, Is.EqualTo(RssSnapshot.Unavailable));
        Assert.That((s.RssAnonKb, s.RssFileKb, s.VmHwmKb, s.Threads), Is.EqualTo((0L, 0L, 0L, 0L)));
        Assert.That(s.RssAnonMb, Is.Zero);
    }

    [TestCase("Name:\tx\nVmRSS:\t 100 kB\n", TestName = "a kernel before 4.5: no RssAnon line")]
    [TestCase("", TestName = "empty")]
    [TestCase("RssAnon: kB\n", TestName = "RssAnon without a number")]
    [TestCase("not a status file at all", TestName = "garbage")]
    public void WithoutAnAnonLine_IsUnavailable(string text)
    {
        Assert.That(new RssReader(new FakeProcSource().With("/proc/self/status", text)).Read().Available, Is.False);
    }

    [Test]
    public void AMissingOptionalLine_ReadsAsZero()
    {
        RssSnapshot s = new RssReader(new FakeProcSource().With("/proc/self/status", "RssAnon: 10 kB\n")).Read();
        Assert.That(s, Is.EqualTo(new RssSnapshot(true, 10, 0, 0, 0)));
    }

    [Test]
    public void AHostileSource_IsUnavailable_NotAnException()
    {
        var source = new FakeProcSource().With("/proc/self/status", Status);
        source.Throws = new IOException("disk on fire");
        Assert.That(new RssReader(source).Read().Available, Is.False);
        source.Throws = new UnauthorizedAccessException("denied");
        Assert.That(new RssReader(source).Read().Available, Is.False);
    }

    [Test]
    public void TheRealSource_WorksOnLinux_AndIsUnavailableElsewhereWithoutThrowing()
    {
        RssSnapshot s = new RssReader().Read();
        if (File.Exists("/proc/self/status") && File.ReadAllText("/proc/self/status").Contains("RssAnon:", StringComparison.Ordinal))
        {
            Assert.That(s.Available, Is.True);
            Assert.That(s.RssAnonKb, Is.GreaterThan(1000), "a .NET process holds more than 1 MB of anonymous memory");
            Assert.That(s.Threads, Is.GreaterThanOrEqualTo(1));
            Assert.That(s.VmHwmKb, Is.GreaterThanOrEqualTo(s.RssAnonKb));
        }
        else
        {
            Assert.That(s.Available, Is.False);
        }
    }

    [Test]
    public void ASourceWithAShortBuffer_StillGetsTheFirstLines()
    {
        // the status file can be long on a busy kernel: the reader's buffer is large, but a source may deliver less; the fields it saw are read
        string longStatus = "RssAnon: 55 kB\n" + new string('x', 20000);
        Assert.That(new RssReader(new FakeProcSource().With("/proc/self/status", longStatus)).Read(), Is.EqualTo(new RssSnapshot(true, 55, 0, 0, 0)));
    }
}

[TestFixture]
public class RssSamplerTests
{
    [Test]
    public void KeepsTheMaxima_OfEverySample_AndTheLast()
    {
        // every read of the fake serves the next value: 100, 900, 300, 50, ... (the sampler reads once at Start, once per interval, once at Stop)
        int[] anon = [100, 900, 300, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50];
        var source = new FakeProcSource().WithSeries("/proc/self/status", call => $"RssAnon: {anon[Math.Min(call, anon.Length - 1)]} kB\nRssFile: {call} kB\nVmHWM: {1000 + call} kB\nThreads: {1 + (call % 5)}\n");
        using var sampler = new RssSampler(new RssReader(source), 1);
        sampler.Start();
        Assert.That(SpinWait.SpinUntil(() => source.Calls >= 12, 10_000), Is.True, "the sampler keeps polling");
        RssStats stats = sampler.Stop();
        Assert.That(stats.Available, Is.True);
        Assert.That(stats.Samples, Is.GreaterThanOrEqualTo(3));
        Assert.That(stats.MaxRssAnonKb, Is.EqualTo(900), "the maximum, not the last or the first");
        Assert.That(stats.LastRssAnonKb, Is.EqualTo(50), "the last sample (taken at Stop)");
        Assert.That(stats.MaxRssFileKb, Is.GreaterThanOrEqualTo(2));
        Assert.That(stats.MaxVmHwmKb, Is.GreaterThanOrEqualTo(1002));
        Assert.That(stats.MaxThreads, Is.GreaterThanOrEqualTo(3));
    }

    [Test]
    public void Unavailable_StartsNoThread_ReportsUnavailable_AndNeverThrows()
    {
        var source = new FakeProcSource();
        using var sampler = new RssSampler(new RssReader(source), 1);
        sampler.Start();
        RssStats stats = sampler.Stop();
        Assert.That(stats.Available, Is.False);
        Assert.That(stats.Samples, Is.Zero);
        Assert.That(stats.MaxRssAnonKb, Is.Zero);
        Assert.That(source.Calls, Is.LessThanOrEqualTo(3), "no polling loop where there is nothing to poll");
        Assert.That(sampler.Stop().Available, Is.False, "a second Stop is harmless");
    }

    [Test]
    public void StopWithoutStart_ReturnsEmptyStats_StartTwiceIsRefused_DisposeStops()
    {
        var source = new FakeProcSource().With("/proc/self/status", "RssAnon: 5 kB\n");
        var idle = new RssSampler(new RssReader(source), 5);
        RssStats none = idle.Stop();
        Assert.That(none.Samples, Is.Zero);
        idle.Dispose();
        var s = new RssSampler(new RssReader(source), 5);
        s.Start();
        Assert.Throws<InvalidOperationException>(s.Start);
        s.Dispose(); // joins the thread
        s.Dispose();
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new RssSampler(new RssReader(source), 0));
        Assert.Throws<ArgumentNullException>(() => _ = new RssSampler(null!, 5));
    }

    [Test]
    public void TheRealSampler_SeesAnAllocation()
    {
        if (!new RssReader().Read().Available) Assert.Ignore("no /proc/self/status here");
        using var sampler = new RssSampler(new RssReader(), 2);
        sampler.Start();
        long before = new RssReader().Read().RssAnonKb;
        var block = new byte[64 * 1024 * 1024];
        for (int i = 0; i < block.Length; i += 4096) block[i] = 1;
        RssStats stats = sampler.Stop(); // the last sample is taken at Stop, while the block is alive
        GC.KeepAlive(block);
        Assert.That(stats.MaxRssAnonKb, Is.GreaterThanOrEqualTo(before + (50 * 1024)), "a touched 64 MB array shows in RssAnon");
    }
}
