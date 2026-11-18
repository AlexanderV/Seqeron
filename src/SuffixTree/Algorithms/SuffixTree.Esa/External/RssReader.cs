using System.Diagnostics;

namespace SuffixTree.Esa.External;

/// <summary>One reading of <c>/proc/self/status</c> in KiB as the kernel prints it; <see cref="Available"/> is false where the file or its <c>RssAnon</c> line does not exist.</summary>
internal readonly record struct RssSnapshot(bool Available, long RssAnonKb, long RssFileKb, long VmHwmKb, long Threads)
{
    /// <summary>The reading of a platform without the information.</summary>
    internal static RssSnapshot Unavailable => default;

    /// <summary>Anonymous resident memory in MB (10^6 bytes).</summary>
    internal double RssAnonMb => RssAnonKb * 1024.0 / 1e6;
}

/// <summary>
/// Reads <c>RssAnon</c> (anonymous resident pages: the quantity the working-set cap of the builder is gated on, ADR-004 5 and the owner decision of 2026-11-18), <c>RssFile</c> (mapped file pages, reported only),
/// <c>VmHWM</c> and <c>Threads</c> from <c>/proc/self/status</c>. Never throws: where the file or the <c>RssAnon</c> line is missing the reading is <see cref="RssSnapshot.Unavailable"/>.
/// </summary>
internal sealed class RssReader
{
    private const string Path = "/proc/self/status";
    private readonly IProcSource _source;

    internal RssReader(IProcSource? source = null) => _source = source ?? FileProcSource.Instance;

    internal RssSnapshot Read()
    {
        Span<byte> buffer = stackalloc byte[8192];
        int got = ProcText.SafeRead(_source, Path, buffer);
        if (got <= 0) return RssSnapshot.Unavailable;
        ReadOnlySpan<byte> text = buffer[..got];
        if (!ProcText.TryGetLong(text, "RssAnon:"u8, out long anon)) return RssSnapshot.Unavailable;
        ProcText.TryGetLong(text, "RssFile:"u8, out long file);
        ProcText.TryGetLong(text, "VmHWM:"u8, out long hwm);
        ProcText.TryGetLong(text, "Threads:"u8, out long threads);
        return new RssSnapshot(true, anon, file, hwm, threads);
    }
}

/// <summary>Maxima observed by an <see cref="RssSampler"/>, in KiB as printed by the kernel.</summary>
internal readonly record struct RssStats(bool Available, long Samples, long MaxRssAnonKb, long MaxRssFileKb, long MaxVmHwmKb, long MaxThreads, long LastRssAnonKb);

/// <summary>
/// Background sampler of an <see cref="RssReader"/>: every <c>intervalMs</c> it reads the status and keeps the maxima. One immediate sample is taken by <see cref="Start"/> and a last one by <see cref="Stop"/>, so a short
/// run still has data. Where the reader is unavailable no thread is started, <see cref="Stop"/> reports <see cref="RssStats.Available"/> = false and nothing throws.
/// </summary>
internal sealed class RssSampler : IDisposable
{
    private readonly RssReader _reader;
    private readonly int _intervalMs;
    private readonly ManualResetEventSlim _stopSignal = new(false);
    private Thread? _thread;
    private volatile bool _stop;
    private bool _started;
    private bool _available;
    private long _samples;
    private long _maxAnon;
    private long _maxFile;
    private long _maxHwm;
    private long _maxThreads;
    private long _lastAnon;

    internal RssSampler(RssReader reader, int intervalMs = 10)
    {
        ArgumentNullException.ThrowIfNull(reader);
        ArgumentOutOfRangeException.ThrowIfLessThan(intervalMs, 1);
        _reader = reader;
        _intervalMs = intervalMs;
    }

    /// <summary>Takes the first sample and, when the platform provides one, starts the sampler thread.</summary>
    internal void Start()
    {
        if (_started) throw new InvalidOperationException("sampler already started");
        _started = true;
        _stop = false;
        _stopSignal.Reset();
        RssSnapshot first = _reader.Read();
        if (!first.Available) return;
        _available = true;
        Record(first);
        _thread = new Thread(Loop) { IsBackground = true, Name = "RssSampler", Priority = ThreadPriority.AboveNormal };
        _thread.Start();
    }

    /// <summary>Stops the thread, takes a last sample and returns the maxima; empty statistics when never started.</summary>
    internal RssStats Stop()
    {
        if (_thread is not null)
        {
            _stop = true;
            _stopSignal.Set(); // wakes the thread from its wait at once, so a long interval does not delay the stop
            _thread.Join();
            _thread = null;
            RssSnapshot last = _reader.Read();
            if (last.Available) Record(last);
        }
        _started = false;
        return new RssStats(_available, Interlocked.Read(ref _samples), Interlocked.Read(ref _maxAnon), Interlocked.Read(ref _maxFile), Interlocked.Read(ref _maxHwm), Interlocked.Read(ref _maxThreads), Interlocked.Read(ref _lastAnon));
    }

    public void Dispose()
    {
        if (_started) Stop();
        _stopSignal.Dispose();
    }

    private void Loop()
    {
        long start = Stopwatch.GetTimestamp();
        long period = Stopwatch.Frequency * _intervalMs / 1000;
        long next = start + period;
        while (!_stop)
        {
            long now = Stopwatch.GetTimestamp();
            if (now < next)
            {
                int sleep = (int)((next - now) * 1000 / Stopwatch.Frequency);
                if (sleep > 0) _stopSignal.Wait(sleep);
                else Thread.Yield();
                continue;
            }
            RssSnapshot s = _reader.Read();
            if (s.Available) Record(s);
            next += period;
            if (next < now) next = now + period; // fell behind (a long pause): do not burst
        }
    }

    private void Record(RssSnapshot s)
    {
        Interlocked.Increment(ref _samples);
        Interlocked.Exchange(ref _lastAnon, s.RssAnonKb);
        Atomic.RaiseTo(ref _maxAnon, s.RssAnonKb);
        Atomic.RaiseTo(ref _maxFile, s.RssFileKb);
        Atomic.RaiseTo(ref _maxHwm, s.VmHwmKb);
        Atomic.RaiseTo(ref _maxThreads, s.Threads);
    }
}
