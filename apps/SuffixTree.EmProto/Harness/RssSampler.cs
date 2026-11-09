using System.Diagnostics;
using System.Text;
using Microsoft.Win32.SafeHandles;

namespace SuffixTree.EmProto.Harness;

/// <summary>Maxima observed by an <see cref="RssSampler"/>, in KiB as printed by the kernel (E0 of ADR-004).</summary>
internal readonly record struct RssStats(long Samples, long MaxRssAnonKb, long MaxRssFileKb, long MaxVmHwmKb, long MaxThreads, long LastRssAnonKb)
{
    /// <summary>Peak anonymous resident memory in MB (10^6 bytes).</summary>
    internal double MaxRssAnonMb => MaxRssAnonKb * 1024.0 / 1e6;

    /// <summary>Peak resident file-backed (mapped) memory in MB.</summary>
    internal double MaxRssFileMb => MaxRssFileKb * 1024.0 / 1e6;
}

/// <summary>
/// Background sampler of <c>/proc/self/status</c> (ADR-004 section 5, E0): every <c>interval</c> ms it reads <c>RssAnon</c> (anonymous resident pages, the quantity
/// that the 256 MB cap of assumption A1 is gated on), <c>RssFile</c> (mapped file pages, reported only), <c>VmHWM</c> and <c>Threads</c>, and keeps the maxima.
/// One dedicated thread, one open handle, one reused buffer: no allocation per sample. Linux only (<see cref="IsSupported"/>).
/// </summary>
internal sealed class RssSampler : IDisposable
{
    private static readonly byte[] KeyRssAnon = Encoding.ASCII.GetBytes("RssAnon:");
    private static readonly byte[] KeyRssFile = Encoding.ASCII.GetBytes("RssFile:");
    private static readonly byte[] KeyHwm = Encoding.ASCII.GetBytes("VmHWM:");
    private static readonly byte[] KeyThreads = Encoding.ASCII.GetBytes("Threads:");

    private readonly int _intervalMs;
    private readonly string _path;
    private readonly List<(double Seconds, long RssAnonKb)>? _series;
    private Thread? _thread;
    private volatile bool _stop;
    private long _samples;
    private long _maxAnon;
    private long _maxFile;
    private long _maxHwm;
    private long _maxThreads;
    private long _lastAnon;

    internal RssSampler(int intervalMs = 10, string path = "/proc/self/status", bool keepSeries = false)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(intervalMs, 1);
        _intervalMs = intervalMs;
        _path = path;
        if (keepSeries) _series = [];
    }

    /// <summary>True when <c>/proc/self/status</c> exists and has an <c>RssAnon</c> line (Linux 4.5 and later).</summary>
    internal static bool IsSupported
    {
        get
        {
            try
            {
                return File.Exists("/proc/self/status") && File.ReadAllText("/proc/self/status").Contains("RssAnon:", StringComparison.Ordinal);
            }
            catch (IOException)
            {
                return false;
            }
        }
    }

    /// <summary>The (seconds since start, RssAnon KiB) series when requested at construction, else empty.</summary>
    internal IReadOnlyList<(double Seconds, long RssAnonKb)> Series => _series ?? (IReadOnlyList<(double, long)>)[];

    /// <summary>Starts the sampler thread (one immediate sample is taken before it returns, so a short run still has data).</summary>
    internal void Start()
    {
        if (_thread is not null) throw new InvalidOperationException("sampler already started");
        _stop = false;
        using SafeFileHandle first = File.OpenHandle(_path);
        SampleOnce(first, new byte[8192], 0.0);
        _thread = new Thread(Loop) { IsBackground = true, Name = "RssSampler", Priority = ThreadPriority.AboveNormal };
        _thread.Start();
    }

    /// <summary>Stops the thread, takes a last sample and returns the maxima.</summary>
    internal RssStats Stop()
    {
        if (_thread is null) return Snapshot();
        _stop = true;
        _thread.Join();
        _thread = null;
        using SafeFileHandle last = File.OpenHandle(_path);
        SampleOnce(last, new byte[8192], double.NaN);
        return Snapshot();
    }

    public void Dispose()
    {
        if (_thread is not null) Stop();
    }

    private RssStats Snapshot()
        => new(Interlocked.Read(ref _samples), Interlocked.Read(ref _maxAnon), Interlocked.Read(ref _maxFile), Interlocked.Read(ref _maxHwm), Interlocked.Read(ref _maxThreads), Interlocked.Read(ref _lastAnon));

    private void Loop()
    {
        using SafeFileHandle handle = File.OpenHandle(_path);
        var buffer = new byte[8192];
        long start = Stopwatch.GetTimestamp();
        long period = Stopwatch.Frequency * _intervalMs / 1000;
        long next = start + period;
        while (!_stop)
        {
            long now = Stopwatch.GetTimestamp();
            if (now < next)
            {
                int sleep = (int)((next - now) * 1000 / Stopwatch.Frequency);
                if (sleep > 0) Thread.Sleep(sleep);
                else Thread.Yield();
                continue;
            }
            SampleOnce(handle, buffer, (now - start) / (double)Stopwatch.Frequency);
            next += period;
            if (next < now) next = now + period; // fell behind (long pause): do not burst
        }
    }

    private void SampleOnce(SafeFileHandle handle, byte[] buffer, double seconds)
    {
        int length = RandomAccess.Read(handle, buffer, 0);
        ReadOnlySpan<byte> text = buffer.AsSpan(0, length);
        long anon = Field(text, KeyRssAnon);
        long file = Field(text, KeyRssFile);
        long hwm = Field(text, KeyHwm);
        long threads = Field(text, KeyThreads);
        Interlocked.Increment(ref _samples);
        Interlocked.Exchange(ref _lastAnon, anon);
        Atomic.RaiseTo(ref _maxAnon, anon);
        Atomic.RaiseTo(ref _maxFile, file);
        Atomic.RaiseTo(ref _maxHwm, hwm);
        Atomic.RaiseTo(ref _maxThreads, threads);
        if (_series is not null && !double.IsNaN(seconds))
        {
            lock (_series) _series.Add((seconds, anon));
        }
    }

    private static long Field(ReadOnlySpan<byte> text, ReadOnlySpan<byte> key)
    {
        int at = text.IndexOf(key);
        if (at < 0) return 0;
        ReadOnlySpan<byte> rest = text[(at + key.Length)..];
        int i = 0;
        while (i < rest.Length && (rest[i] == (byte)' ' || rest[i] == (byte)'\t')) i++;
        long value = 0;
        while (i < rest.Length && rest[i] >= (byte)'0' && rest[i] <= (byte)'9')
        {
            value = (value * 10) + (rest[i] - (byte)'0');
            i++;
        }
        return value;
    }
}
