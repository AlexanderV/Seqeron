using System.Globalization;

namespace SuffixTree.EmProto.Harness;

/// <summary>Maxima of the system-wide dirty page-cache bytes and of a cgroup's memory usage seen by a <see cref="MemSampler"/>.</summary>
internal readonly record struct MemStats(long MaxDirtyKb, long MaxCgroupUsageBytes, long MaxCachedKb, long CgSamples, long MaxCgRssBytes, double AvgCgRssBytes, double AvgCgMappedFileBytes, double AvgCgCacheBytes);

/// <summary>
/// Samples <c>/proc/meminfo</c> (<c>Dirty</c>, <c>Cached</c>) and, when given, a cgroup v1 <c>memory.usage_in_bytes</c> every 50 ms (E3: dirty pages of the mapped output windows, page-cache use
/// under a limit). Reported, not gated.
/// </summary>
internal sealed class MemSampler(string? cgroupUsagePath = null) : IDisposable
{
    private volatile bool _stop;
    private Thread? _thread;
    private long _maxDirty;
    private long _maxCg;
    private long _maxCached;
    private long _cgSamples;
    private long _maxCgRss;
    private double _sumRss;
    private double _sumMapped;
    private double _sumCache;
    private readonly object _gate = new();

    internal void Start()
    {
        Sample();
        _thread = new Thread(() =>
        {
            while (!_stop)
            {
                Sample();
                Thread.Sleep(50);
            }
        })
        { IsBackground = true, Name = "MemSampler" };
        _thread.Start();
    }

    internal MemStats Stop()
    {
        _stop = true;
        _thread?.Join();
        _thread = null;
        Sample();
        lock (_gate)
        {
            long k = Math.Max(1, _cgSamples);
            return new MemStats(Interlocked.Read(ref _maxDirty), Interlocked.Read(ref _maxCg), Interlocked.Read(ref _maxCached), _cgSamples, Interlocked.Read(ref _maxCgRss), _sumRss / k, _sumMapped / k, _sumCache / k);
        }
    }

    public void Dispose()
    {
        if (_thread is not null) Stop();
    }

    private void Sample()
    {
        try
        {
            foreach (string line in File.ReadLines("/proc/meminfo"))
            {
                if (line.StartsWith("Dirty:", StringComparison.Ordinal)) Atomic.RaiseTo(ref _maxDirty, Kb(line));
                else if (line.StartsWith("Cached:", StringComparison.Ordinal)) Atomic.RaiseTo(ref _maxCached, Kb(line));
            }
            if (cgroupUsagePath is not null)
            {
                Atomic.RaiseTo(ref _maxCg, long.Parse(File.ReadAllText(cgroupUsagePath).Trim(), CultureInfo.InvariantCulture));
                string stat = Path.Combine(Path.GetDirectoryName(cgroupUsagePath)!, "memory.stat");
                long rss = 0, mapped = 0, cache = 0;
                foreach (string line in File.ReadLines(stat))
                {
                    string[] parts = line.Split(' ');
                    if (parts.Length != 2) continue;
                    if (parts[0] == "rss") rss = long.Parse(parts[1], CultureInfo.InvariantCulture);
                    else if (parts[0] == "mapped_file") mapped = long.Parse(parts[1], CultureInfo.InvariantCulture);
                    else if (parts[0] == "cache") cache = long.Parse(parts[1], CultureInfo.InvariantCulture);
                }
                Atomic.RaiseTo(ref _maxCgRss, rss);
                lock (_gate)
                {
                    _cgSamples++;
                    _sumRss += rss;
                    _sumMapped += mapped;
                    _sumCache += cache;
                }
            }
        }
        catch (IOException)
        {
            // a sample may fail while the cgroup is being torn down; the next one will do
        }
    }

    private static long Kb(string line)
    {
        string[] parts = line.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries);
        return long.Parse(parts[1], CultureInfo.InvariantCulture);
    }
}
