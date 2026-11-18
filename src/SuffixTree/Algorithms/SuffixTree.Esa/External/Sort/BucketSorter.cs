using System.Buffers;
using System.Diagnostics;
using System.Runtime.ExceptionServices;

namespace SuffixTree.Esa.External;

/// <summary>Parameters of PH2 and PH3 (ADR-004 sections 4.2 and 5); the tests scale them down to K = 3..50, F = 2..4, X = 7..64 so that every multi-level partition and every spill path runs on texts of 10^1-10^4 symbols. The output never depends on any of them.</summary>
internal sealed record SortOptions
{
    /// <summary>K: records per leaf bucket (a bucket of at most K positions is sorted in RAM).</summary>
    internal long LeafCapacity { get; init; } = 1 << 20;

    /// <summary>F: largest number of children of one partition step (the fan-out of ADR-004: 64).</summary>
    internal int FanOut { get; init; } = 64;

    /// <summary>Random sample positions per child when splitters are chosen.</summary>
    internal int Oversample { get; init; } = 16;

    /// <summary>Worker threads (the calling thread is one of them).</summary>
    internal int Threads { get; init; } = 1;

    /// <summary>Seed of the splitter sampling; the output does not depend on it, the shape of the partition does.</summary>
    internal ulong Seed { get; init; } = 1;

    /// <summary>Directory of the position files of the buckets (created when missing; every file of a sort is removed when it ends).</summary>
    internal required string ScratchDir { get; init; }

    /// <summary>Write buffer per child bucket (a bucket that fits stays in memory).</summary>
    internal int SpillBytes { get; init; } = 64 * 1024;

    /// <summary>Prefix of the names of the position files (a sort that runs inside another build gives itself its own prefix, so that the failure path of the outer build finds its files).</summary>
    internal string FilePrefix { get; init; } = "pos";

    /// <summary>The ledger of the scratch bytes (null = a private one over <see cref="Counters"/>).</summary>
    internal ScratchLedger? Ledger { get; init; }

    /// <summary>Shared cap of the native working buffers (null = unlimited).</summary>
    internal ArenaBudget? Budget { get; init; }

    /// <summary>Test hook: options of the arenas (guard bands, poison pattern); null = production arenas.</summary>
    internal ArenaOptions? ArenaGuards { get; init; }

    /// <summary>Receives the explicit I/O and the scratch bytes (null = not counted).</summary>
    internal IoCounters? Counters { get; init; }

    /// <summary>Test hook, called with <c>"partition"</c> before a bucket is split and <c>"leaf"</c> before a leaf is sorted (from the worker's thread); an exception thrown here fails the sort like a failing step.</summary>
    internal Action<string>? Trace { get; init; }
}

/// <summary>What a sort did, for the report.</summary>
internal sealed class SortStats
{
    internal SortWork Work { get; } = new();

    internal long Buckets { get; set; }

    internal long Leaves { get; set; }

    internal long Partitions { get; set; }

    internal int MaxDepth { get; set; }

    internal long MaxLeaf { get; set; }

    /// <summary>Degenerate splits (all positions on one side) that were repeated with more samples.</summary>
    internal long Retries { get; set; }

    internal double PartitionSeconds { get; set; }

    internal double LeafSeconds { get; set; }
}

/// <summary>One bucket still to be processed: its positions, the entry index of its first suffix in the final suffix array, and a deterministic id (seeds the sampling).</summary>
internal sealed record BuildTask(IPosSource Bucket, long SaOffset, ulong Id, int Depth);

/// <summary>
/// The bucket partition and the bucket sort of ADR-004 (PH2 and PH3, with the comparator of the difference-cover sample): the suffixes are distributed over lexicographic buckets by splitters (real suffixes ordered by the
/// comparator, a stratified random sample), every bucket larger than K is split again, every bucket of at most K positions is sorted in RAM and its slice of the suffix array is written at its final offset by the writer of the
/// worker. The output is the unique suffix array: independent of the thread count, K, F and the seed. A failure of any worker ends the sort; the exception reaches the caller with its own type, and the position files of the sort,
/// the writers and the arenas are given back.
/// </summary>
internal sealed class BucketSorter
{
    private readonly MappedTextReader _text;
    private readonly SortOptions _o;
    private readonly ISaSliceSink _sink;
    private readonly RankTable _ranks;
    private readonly ArenaBudget _budget;
    private readonly PositionFiles _files;
    private readonly WorkPool<BuildTask> _pool = new();
    private readonly SortStats _stats = new();
    private readonly object _statsGate = new();

    private BucketSorter(MappedTextReader text, SortOptions options, ISaSliceSink sink, RankTable ranks, PositionFiles files)
    {
        _text = text;
        _o = options;
        _sink = sink;
        _ranks = ranks;
        _budget = options.Budget ?? new ArenaBudget(long.MaxValue);
        _files = files;
    }

    /// <summary>Sorts all suffixes of <paramref name="text"/> and writes the suffix array (<c>n + 1</c> entries, <c>SA[0] = n</c>) through <paramref name="sink"/>.</summary>
    internal static SortStats Run(MappedTextReader text, SortOptions options, ISaSliceSink sink, RankTable ranks)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(options);
        ArgumentNullException.ThrowIfNull(sink);
        ArgumentNullException.ThrowIfNull(ranks);
        ArgumentNullException.ThrowIfNull(options.ScratchDir);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.LeafCapacity, 2);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(options.LeafCapacity, int.MaxValue / 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.FanOut, 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.Threads, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.Oversample, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.SpillBytes, 1);
        if (text.Length > PackedTextFormat.MaxSymbols) throw new NotSupportedException($"{text.Length} symbols: positions are 32-bit, the limit is {PackedTextFormat.MaxSymbols}");
        string dir = Path.GetFullPath(options.ScratchDir);
        bool createdRoot = !Directory.Exists(dir);
        ScratchLedger ledger = options.Ledger ?? new ScratchLedger(options.Counters);
        var files = new PositionFiles(dir, options.FilePrefix, ledger, options.Counters);
        try
        {
            Directory.CreateDirectory(dir);
            return new BucketSorter(text, options, sink, ranks, files).Execute();
        }
        catch
        {
            if (options.Ledger is null) ledger.ReleaseAll();
            DeleteFilesQuietly(dir, options.FilePrefix + "-*");
            throw;
        }
        finally
        {
            if (createdRoot) DeleteEmptyDirectoryQuietly(dir);
        }
    }

    private SortStats Execute()
    {
        long n = _text.Length;
        using ISaSliceWriter main = _sink.CreateWriter();
        main.Write(0, [(uint)n]);
        if (n == 0) return _stats;
        _pool.Push(new BuildTask(new RangePositions(n), 1, 1UL, 0));
        int threads = Math.Max(1, _o.Threads);
        var extra = new List<Thread>();
        for (int t = 1; t < threads; t++)
        {
            var th = new Thread(() => Worker(null)) { IsBackground = true, Name = $"ph3-worker-{t}" };
            extra.Add(th);
            th.Start();
        }
        Worker(main);
        foreach (Thread th in extra) th.Join();
        if (_pool.Error is { } error) ExceptionDispatchInfo.Capture(error).Throw();
        return _stats;
    }

    /// <summary>One worker: the calling thread brings its writer, every other thread makes its own (and disposes it).</summary>
    private void Worker(ISaSliceWriter? given)
    {
        ISaSliceWriter? own = null;
        try
        {
            if (given is null) own = _sink.CreateWriter();
            RunTasks(given ?? own!);
        }
#pragma warning disable CA1031 // a failure of any kind in a worker thread is recorded and rethrown, with its own type and stack, on the calling thread (an unhandled exception would end the process)
        catch (Exception e)
        {
            _pool.Fail(e);
        }
#pragma warning restore CA1031
        finally
        {
            own?.Dispose();
        }
    }

    private void RunTasks(ISaSliceWriter writer)
    {
        var work = new SortWork();
        using var arena = new BudgetedArena(_budget, _o.ArenaGuards);
        var order = new SuffixOrder(_text, _ranks, work);
        LeafSorter? leaf = null;
        long buckets = 0, leaves = 0, partitions = 0, retries = 0, maxLeaf = 0;
        int maxDepth = 0;
        long partTicks = 0, leafTicks = 0;
        while (_pool.Take(out BuildTask task))
        {
            try
            {
                buckets++;
                maxDepth = Math.Max(maxDepth, task.Depth);
                long t0 = Stopwatch.GetTimestamp();
                if (task.Bucket.Count <= _o.LeafCapacity)
                {
                    _o.Trace?.Invoke("leaf");
                    leaf ??= new LeafSorter(arena, (int)_o.LeafCapacity);
                    leaf.Sort(task.Bucket, task.SaOffset, order, writer);
                    task.Bucket.Dispose();
                    leaves++;
                    maxLeaf = Math.Max(maxLeaf, task.Bucket.Count);
                    leafTicks += Stopwatch.GetTimestamp() - t0;
                }
                else
                {
                    _o.Trace?.Invoke("partition");
                    retries += Partition(task, order, work);
                    partitions++;
                    partTicks += Stopwatch.GetTimestamp() - t0;
                }
            }
#pragma warning disable CA1031 // see Worker
            catch (Exception e)
            {
                task.Bucket.Dispose(); // the file of the bucket and its descriptor are given back
                _pool.Fail(e);
            }
#pragma warning restore CA1031
            finally
            {
                _pool.Done();
            }
        }
        arena.VerifyCanaries();
        lock (_statsGate)
        {
            _stats.Work.Add(work);
            _stats.Buckets += buckets;
            _stats.Leaves += leaves;
            _stats.Partitions += partitions;
            _stats.Retries += retries;
            _stats.MaxDepth = Math.Max(_stats.MaxDepth, maxDepth);
            _stats.MaxLeaf = Math.Max(_stats.MaxLeaf, maxLeaf);
            _stats.PartitionSeconds += partTicks / (double)Stopwatch.Frequency;
            _stats.LeafSeconds += leafTicks / (double)Stopwatch.Frequency;
        }
    }

    /// <summary>Splits an oversized bucket into children by splitters; returns the number of degenerate attempts that had to be repeated.</summary>
    private long Partition(BuildTask task, SuffixOrder order, SortWork work)
    {
        IPosSource src = task.Bucket;
        long m = src.Count;
        long retries = 0;
        int oversample = _o.Oversample;
        for (int attempt = 0; attempt < 64; attempt++)
        {
            var rng = new SeededRandom(_o.Seed ^ (task.Id * 0x9E3779B97F4A7C15UL) ^ ((ulong)attempt * 0xD1B54A32D192ED03UL));
            int fan = (int)Math.Min(Math.Min(_o.FanOut, m), Math.Max(2, ((2 * m) + _o.LeafCapacity - 1) / _o.LeafCapacity));
            int sampleCount = (int)Math.Min(m, (long)fan * oversample);
            var sample = new uint[sampleCount];
            for (int t = 0; t < sampleCount; t++)
            {
                long lo = t * m / sampleCount;
                long hi = (t + 1) * m / sampleCount;
                sample[t] = src.At(lo + rng.NextInt((int)(hi - lo)));
            }
            sample.AsSpan().Sort(new FullComparer(order));
            int nSpl = fan - 1;
            var splPos = new uint[nSpl];
            var splKey = new ulong[nSpl];
            for (int t = 1; t < fan; t++)
            {
                splPos[t - 1] = sample[(long)t * sampleCount / fan];
                splKey[t - 1] = _text.Key(splPos[t - 1]);
            }

            var writers = new PosWriter[fan];
            uint[]? chunk = null;
            bool handedOver = false;
            try
            {
                for (int b = 0; b < fan; b++) writers[b] = new PosWriter(_files, _o.SpillBytes, (int)Math.Min((m / fan) + 16, 1 << 16));
                chunk = ArrayPool<uint>.Shared.Rent(1 << 14);
                for (long start = 0; start < m; start += 1 << 14)
                {
                    int take = (int)Math.Min(1 << 14, m - start);
                    src.Read(start, chunk.AsSpan(0, take));
                    for (int t = 0; t < take; t++)
                    {
                        uint x = chunk[t];
                        ulong kx = _text.Key(x);
                        int lo = 0, hi = nSpl;
                        while (lo < hi)
                        {
                            int mid = (lo + hi) >> 1;
                            ulong ks = splKey[mid];
                            if (ks < kx || (ks == kx && order.CompareTied(splPos[mid], x) <= 0)) lo = mid + 1;
                            else hi = mid;
                        }
                        writers[lo].Append(x);
                    }
                }
                work.ScanPositions += m;
                long biggest = 0;
                foreach (PosWriter w in writers) biggest = Math.Max(biggest, w.Count);
                if (biggest == m)
                {
                    // degenerate split (all positions on one side): sample more and try again; with sampleCount == m the splitters are exact quantiles, so this ends
                    retries++;
                    oversample *= 2;
                    continue;
                }
                long offset = task.SaOffset;
                var children = new List<BuildTask>(fan);
                try
                {
                    for (int b = 0; b < fan; b++)
                    {
                        long count = writers[b].Count;
                        IPosSource child = writers[b].Complete();
                        if (count == 0)
                        {
                            child.Dispose();
                            continue;
                        }
                        children.Add(new BuildTask(child, offset, (task.Id * 0x100000001B3UL) ^ (ulong)(b + 1), task.Depth + 1));
                        offset += count;
                    }
                    handedOver = true;
                }
                finally
                {
                    if (!handedOver) foreach (BuildTask child in children) child.Bucket.Dispose();
                }
                src.Dispose();
                for (int b = children.Count - 1; b >= 0; b--) _pool.Push(children[b]);
                return retries;
            }
            finally
            {
                foreach (PosWriter? w in writers) w?.Dispose();
                if (chunk is not null) ArrayPool<uint>.Shared.Return(chunk);
            }
        }
        throw new InvalidOperationException("the partition of a bucket did not shrink in 64 attempts");
    }

    private static void DeleteFilesQuietly(string dir, string pattern)
    {
        try
        {
            if (!Directory.Exists(dir)) return;
            foreach (string file in Directory.GetFiles(dir, pattern)) RunStore.DeleteQuietly(file);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // reported by the files that are left
        }
    }

    private static void DeleteEmptyDirectoryQuietly(string path)
    {
        try
        {
            if (Directory.Exists(path) && !Directory.EnumerateFileSystemEntries(path).Any()) Directory.Delete(path);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // someone else put a file there between the check and the delete: it is theirs
        }
    }
}
