using System.Buffers;
using System.Diagnostics;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Build;

/// <summary>One bucket still to be processed: its positions, the entry index of its first suffix in the final suffix array, and a deterministic id (seeds the sampling).</summary>
internal sealed record BuildTask(IPosSource Bucket, long SaOffset, ulong Id, int Depth);

/// <summary>LIFO work list with termination detection (depth-first order keeps the number of pending buckets, hence scratch files, small).</summary>
internal sealed class TaskPool
{
    private readonly Stack<BuildTask> _stack = new();
    private int _active;
    private Exception? _error;

    internal Exception? Error => _error;

    internal void Push(BuildTask t)
    {
        lock (_stack)
        {
            _stack.Push(t);
            Monitor.Pulse(_stack);
        }
    }

    /// <summary>Blocks until a task is available (returns true; the caller must call <see cref="Done"/>) or all work is finished or failed (false).</summary>
    internal bool Take(out BuildTask task)
    {
        lock (_stack)
        {
            while (true)
            {
                if (_error is null && _stack.Count > 0)
                {
                    task = _stack.Pop();
                    _active++;
                    return true;
                }
                if (_error is not null || _active == 0)
                {
                    task = null!;
                    Monitor.PulseAll(_stack);
                    return false;
                }
                Monitor.Wait(_stack);
            }
        }
    }

    internal void Done()
    {
        lock (_stack)
        {
            _active--;
            Monitor.PulseAll(_stack);
        }
    }

    internal void Fail(Exception e)
    {
        lock (_stack)
        {
            _error ??= e;
            Monitor.PulseAll(_stack);
        }
    }
}

/// <summary>
/// The bucket partition and the bucket sort of ADR-004 (PH2-PH3, with the comparator of E1 or E2): the suffixes are distributed over lexicographic buckets by splitters (real suffixes
/// ordered by the comparator, a stratified random sample), every bucket larger than K is split again, every bucket of at most K positions is sorted in RAM and its slice of the suffix array is
/// written at its final offset. The output is the unique suffix array: independent of the thread count, K, F and the seed.
/// </summary>
internal sealed class EmBuilder
{
    private readonly PackedText _text;
    private readonly BuildOptions _o;
    private readonly ISaSink _sink;
    private readonly DcsRankTable? _dcs;
    private readonly TaskPool _pool = new();
    private readonly BuildStats _stats = new();
    private readonly object _statsGate = new();

    private EmBuilder(PackedText text, BuildOptions options, ISaSink sink, DcsRankTable? dcs)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(options.LeafCapacity, 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.FanOut, 2);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(options.LeafCapacity, int.MaxValue / 2);
        if (text.Length >= uint.MaxValue) throw new NotSupportedException("positions are 32-bit in the prototype");
        _text = text;
        _o = options;
        _sink = sink;
        _dcs = dcs;
    }

    /// <summary>Builds the suffix array of <paramref name="text"/> (<c>n + 1</c> entries, <c>SA[0] = n</c>) into <paramref name="sink"/>.</summary>
    internal static BuildStats Run(PackedText text, BuildOptions options, ISaSink sink, DcsRankTable? dcs = null)
        => new EmBuilder(text, options, sink, dcs).Execute();

    private BuildStats Execute()
    {
        long n = _text.Length;
        _sink.Write(0, [(uint)n]);
        if (n == 0) return _stats;
        _pool.Push(new BuildTask(new RangePositions(n), 1, 1UL, 0));
        int threads = Math.Max(1, _o.Threads);
        var extra = new List<Thread>();
        for (int t = 1; t < threads; t++)
        {
            var th = new Thread(Worker) { IsBackground = true, Name = $"em-worker-{t}" };
            extra.Add(th);
            th.Start();
        }
        Worker();
        foreach (Thread th in extra) th.Join();
        if (_pool.Error is { } error) throw new InvalidOperationException("build failed: " + error.Message, error);
        return _stats;
    }

    private void Worker()
    {
        var counters = new WorkCounters();
        using var arena = new BudgetedArena(_o.Budget ?? new ArenaBudget(long.MaxValue));
        var order = new SuffixOrder(_text, _dcs, counters);
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
                    leaf ??= new LeafSorter(arena, (int)_o.LeafCapacity, _o.TieSort);
                    leaf.Sort(task.Bucket, task.SaOffset, order, _sink);
                    task.Bucket.Dispose();
                    leaves++;
                    maxLeaf = Math.Max(maxLeaf, task.Bucket.Count);
                    leafTicks += Stopwatch.GetTimestamp() - t0;
                }
                else
                {
                    retries += Partition(task, order, counters);
                    partitions++;
                    partTicks += Stopwatch.GetTimestamp() - t0;
                }
            }
            catch (Exception e)
            {
                _pool.Fail(e);
            }
            finally
            {
                _pool.Done();
            }
        }
        lock (_statsGate)
        {
            _stats.Work.Add(counters);
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
    private long Partition(BuildTask task, SuffixOrder order, WorkCounters counters)
    {
        IPosSource src = task.Bucket;
        long m = src.Count;
        long retries = 0;
        int oversample = _o.Oversample;
        for (int attempt = 0; attempt < 64; attempt++)
        {
            var rng = new SplitMix64(_o.Seed ^ (task.Id * 0x9E3779B97F4A7C15UL) ^ ((ulong)attempt * 0xD1B54A32D192ED03UL));
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
            for (int b = 0; b < fan; b++) writers[b] = new PosWriter(_o.ScratchDir, _o.SpillBytes, (int)Math.Min(m / fan + 16, 1 << 16));
            uint[] chunk = ArrayPool<uint>.Shared.Rent(1 << 14);
            try
            {
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
                counters.ScanPositions += m;
                long biggest = 0;
                foreach (PosWriter w in writers) biggest = Math.Max(biggest, w.Count);
                if (biggest == m)
                {
                    // degenerate split (all positions on one side): sample more and try again; with sampleCount == m the splitters are exact quantiles, so this ends
                    foreach (PosWriter w in writers) w.Dispose();
                    retries++;
                    oversample *= 2;
                    continue;
                }
                long offset = task.SaOffset;
                var children = new List<BuildTask>(fan);
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
                src.Dispose();
                for (int b = children.Count - 1; b >= 0; b--) _pool.Push(children[b]);
                return retries;
            }
            finally
            {
                ArrayPool<uint>.Shared.Return(chunk);
            }
        }
        throw new InvalidOperationException("the partition of a bucket did not shrink in 64 attempts");
    }
}
