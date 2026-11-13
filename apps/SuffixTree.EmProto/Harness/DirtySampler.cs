using System.Diagnostics;
using System.Text;
using Microsoft.Win32.SafeHandles;

namespace SuffixTree.EmProto.Harness;

/// <summary>Maxima of the dirty page-cache bytes seen by a <see cref="DirtySampler"/> (KiB as printed by the kernel).</summary>
internal readonly record struct DirtyStats(long Samples, double BaselineKb, long MaxDirtyKb, double MaxDeltaKb, double MeanDeltaKb, long MaxWritebackKb, long MaxProcessDirtyKb, double MarkBaselineKb, double MaxDeltaAfterMarkKb);

/// <summary>
/// Samples the dirty pages every <c>interval</c> ms (ADR-004 section 5, row E5: "mapped dirty pages &lt;= 64 MB sampled from /proc/meminfo Dirty"): <c>Dirty</c> and <c>Writeback</c> of <c>/proc/meminfo</c> (system-wide) and
/// <c>Private_Dirty + Shared_Dirty</c> of <c>/proc/self/smaps_rollup</c> (the dirty pages this process has mapped: the per-process view). <see cref="Start"/> takes the baseline first (the median of five samples), the delta is
/// <c>Dirty - baseline</c>, so that dirty pages of something else that was already in flight do not count. <see cref="Mark"/> takes a second baseline at a moment of the run (fused mode: the start of the PH5 scan, after
/// PH4 has written its scratch) and keeps the maximum after it separately. One thread, handles opened once, no allocation per sample.
/// </summary>
internal sealed class DirtySampler : IDisposable
{
    private static readonly byte[] KeyDirty = Encoding.ASCII.GetBytes("Dirty:");
    private static readonly byte[] KeyWriteback = Encoding.ASCII.GetBytes("Writeback:");
    private static readonly byte[] KeyPrivate = Encoding.ASCII.GetBytes("Private_Dirty:");
    private static readonly byte[] KeyShared = Encoding.ASCII.GetBytes("Shared_Dirty:");

    private readonly int _intervalMs;
    private readonly string? _filePrefix;
    private readonly List<(double Seconds, long DirtyKb, long ProcessKb)>? _series;
    private Thread? _thread;
    private volatile bool _stop;
    private long _samples;
    private double _baseline;
    private long _maxDirty;
    private long _maxWriteback;
    private long _maxProcess;
    private double _sumDelta;
    private long _markBaseline = -1;
    private long _maxAfterMark = long.MinValue;
    private long _lastDirty;
    private long _phaseDirty;
    private long _phaseProcess;

    /// <param name="intervalMs">Sampling period.</param>
    /// <param name="filePrefix">When given, the per-process view is the sum of <c>Private_Dirty + Shared_Dirty</c> of the mappings in <c>/proc/self/smaps</c> whose path starts with this prefix (the dirty pages of the mapped output files; anonymous memory, which smaps_rollup would add, is left out); null = no per-process view.</param>
    /// <param name="keepSeries">Keep every sample (time, Dirty, per-process dirty).</param>
    internal DirtySampler(int intervalMs = 20, string? filePrefix = null, bool keepSeries = false)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(intervalMs, 1);
        _intervalMs = intervalMs;
        _filePrefix = filePrefix;
        if (keepSeries) _series = [];
    }

    internal IReadOnlyList<(double Seconds, long DirtyKb, long ProcessKb)> Series => _series ?? (IReadOnlyList<(double, long, long)>)[];

    /// <summary>Takes the baseline (five samples, 20 ms apart, median) and starts the sampler thread.</summary>
    internal void Start()
    {
        if (_thread is not null) throw new InvalidOperationException("sampler already started");
        using SafeFileHandle meminfo = File.OpenHandle("/proc/meminfo");
        var buffer = new byte[8192];
        var base5 = new long[5];
        for (int i = 0; i < base5.Length; i++)
        {
            base5[i] = ReadKb(meminfo, buffer, KeyDirty);
            Thread.Sleep(_intervalMs);
        }
        Array.Sort(base5);
        _baseline = base5[2];
        _thread = new Thread(Loop) { IsBackground = true, Name = "DirtySampler", Priority = ThreadPriority.AboveNormal };
        _stop = false;
        _thread.Start();
    }

    /// <summary>The baseline (median of five samples before the start) in KiB; valid after <see cref="Start"/>.</summary>
    internal double BaselineKb => _baseline;

    /// <summary>The maxima of Dirty and of the per-process view since the previous call (or since <see cref="Start"/>), in KiB, and resets them (one call per phase of a pipeline).</summary>
    internal (long DirtyKb, long ProcessKb) TakePhaseMax() => (Interlocked.Exchange(ref _phaseDirty, 0), Interlocked.Exchange(ref _phaseProcess, 0));

    /// <summary>Takes a second baseline now (the current Dirty) and from here on keeps the maximum of the delta to it separately.</summary>
    internal void Mark() => Interlocked.Exchange(ref _markBaseline, Interlocked.Read(ref _lastDirty));

    internal DirtyStats Stop()
    {
        if (_thread is null) return Snapshot();
        _stop = true;
        _thread.Join();
        _thread = null;
        return Snapshot();
    }

    public void Dispose()
    {
        if (_thread is not null) Stop();
    }

    private DirtyStats Snapshot()
    {
        long samples = Math.Max(1, Interlocked.Read(ref _samples));
        long mark = Interlocked.Read(ref _markBaseline);
        long afterMax = Interlocked.Read(ref _maxAfterMark);
        return new DirtyStats(
            samples,
            _baseline,
            Interlocked.Read(ref _maxDirty),
            Interlocked.Read(ref _maxDirty) - _baseline,
            _sumDelta / samples,
            Interlocked.Read(ref _maxWriteback),
            Interlocked.Read(ref _maxProcess),
            mark,
            mark < 0 || afterMax == long.MinValue ? 0 : afterMax - mark);
    }

    private void Loop()
    {
        using SafeFileHandle meminfo = File.OpenHandle("/proc/meminfo");
        using SafeFileHandle? smaps = _filePrefix is not null && File.Exists("/proc/self/smaps") ? File.OpenHandle("/proc/self/smaps") : null;
        var buffer = new byte[8192];
        var smapsBuffer = new byte[1 << 20];
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
            int length = RandomAccess.Read(meminfo, buffer, 0);
            ReadOnlySpan<byte> text = buffer.AsSpan(0, length);
            long dirty = Field(text, KeyDirty);
            long writeback = Field(text, KeyWriteback);
            long process = 0;
            if (smaps is not null) process = FileMappedDirtyKb(smaps, smapsBuffer, _filePrefix!);
            Interlocked.Exchange(ref _lastDirty, dirty);
            Interlocked.Increment(ref _samples);
            Atomic.RaiseTo(ref _maxDirty, dirty);
            Atomic.RaiseTo(ref _maxWriteback, writeback);
            Atomic.RaiseTo(ref _maxProcess, process);
            Atomic.RaiseTo(ref _phaseDirty, dirty);
            Atomic.RaiseTo(ref _phaseProcess, process);
            _sumDelta += dirty - _baseline;
            if (Interlocked.Read(ref _markBaseline) >= 0) Atomic.RaiseTo(ref _maxAfterMark, dirty);
            _series?.Add(((now - start) / (double)Stopwatch.Frequency, dirty, process));
            next += period;
            if (next < now) next = now + period;
        }
    }

    /// <summary>Sum of <c>Private_Dirty</c> and <c>Shared_Dirty</c> of the mappings of <c>/proc/self/smaps</c> whose path starts with <paramref name="prefix"/>.</summary>
    internal static long FileMappedDirtyKb(SafeFileHandle smaps, byte[] buffer, string prefix)
    {
        int total = 0;
        while (total < buffer.Length)
        {
            int got = RandomAccess.Read(smaps, buffer.AsSpan(total), total);
            if (got <= 0) break;
            total += got;
        }
        return SumFileMapped(buffer.AsSpan(0, total), Encoding.UTF8.GetBytes(prefix));
    }

    /// <summary>The parser of <see cref="FileMappedDirtyKb"/>, on the text of an smaps file.</summary>
    internal static long SumFileMapped(ReadOnlySpan<byte> text, ReadOnlySpan<byte> prefix)
    {
        long sum = 0;
        bool inside = false;
        while (!text.IsEmpty)
        {
            int eol = text.IndexOf((byte)'\n');
            ReadOnlySpan<byte> line = eol < 0 ? text : text[..eol];
            text = eol < 0 ? default : text[(eol + 1)..];
            if (line.IsEmpty) continue;
            // a mapping header starts with a hex address range ("7f00-7f10 rw-s 00000000 fe:00 123 /path"); a field line starts with a name and a colon
            int colon = line.IndexOf((byte)':');
            int space = line.IndexOf((byte)' ');
            bool header = space > 0 && (colon < 0 || colon > space || line[..colon].Contains((byte)'-'));
            if (header)
            {
                int slash = line.IndexOf((byte)'/');
                inside = slash >= 0 && line[slash..].StartsWith(prefix);
            }
            else if (inside && (line.StartsWith(KeyPrivate) || line.StartsWith(KeyShared)))
            {
                sum += Field(line, line[..(line.IndexOf((byte)':') + 1)]);
            }
        }
        return sum;
    }

    private static long ReadKb(SafeFileHandle handle, byte[] buffer, byte[] key)
    {
        int length = RandomAccess.Read(handle, buffer, 0);
        return Field(buffer.AsSpan(0, length), key);
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
