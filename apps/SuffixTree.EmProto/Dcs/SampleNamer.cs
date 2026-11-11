using System.Diagnostics;
using SuffixTree.Esa;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Dcs;

/// <summary>
/// The order of the X-prefixes of two sample positions, with the end of the text smaller than every symbol; 0 means equal X-prefixes (both suffixes have at least X symbols and agree on the first X).
/// Exactly the comparison of <see cref="DcsBuilder"/> (so that the names and hence the ranks are the same).
/// </summary>
internal sealed class PrefixOrder(PackedText text, int x)
{
    internal int Compare(long pa, long pb)
    {
        if (pa == pb) return 0;
        long maxLen = text.Length - Math.Max(pa, pb);
        long lim = Math.Min(x, maxLen);
        long steps = 0;
        long l = text.Lce(pa, pb, lim, ref steps);
        if (l >= x) return 0;
        if (l == maxLen) return pa > pb ? -1 : 1;
        return text.Symbol(pa + l) < text.Symbol(pb + l) ? -1 : 1;
    }
}

/// <summary>A sample record used as the lower (inclusive) or upper (exclusive) end of a bucket: the position and its 64-bit key.</summary>
internal sealed record NamingBound(uint Pos, ulong Key);

/// <summary>A node of the bucket tree: a split (children, in order) or a leaf (the result of its sort).</summary>
internal sealed class RangeNode
{
    internal List<RangeNode>? Children { get; set; }

    internal LeafResult? Leaf { get; set; }
}

/// <summary>What a leaf bucket knows after it was sorted and named: counts for the global name offsets and the share of unique names, and its run.</summary>
internal sealed class LeafResult : IDisposable
{
    internal long Count { get; init; }

    /// <summary>Distinct (key, X-prefix) groups in the bucket; the local names are 1..Distinct.</summary>
    internal long Distinct { get; init; }

    /// <summary>Smallest and largest record of the bucket in the order (key, X-prefix, position).</summary>
    internal uint FirstPos { get; init; }

    internal uint LastPos { get; init; }

    /// <summary>Sizes of the first and the last group, and the number of groups of size one.</summary>
    internal long FirstLen { get; init; }

    internal long LastLen { get; init; }

    internal long Singles { get; init; }

    internal BucketRun? Run { get; init; }

    /// <summary>Global name of local name 1, minus 1 (set after all buckets are known).</summary>
    internal uint Offset { get; set; }

    public void Dispose() => Run?.Dispose();
}

/// <summary>The sorted bucket results, in the order of the sample, and the totals.</summary>
internal sealed class NamingResult : IDisposable
{
    internal List<LeafResult> Leaves { get; } = [];

    /// <summary>Distinct names over the whole sample.</summary>
    internal long Distinct { get; set; }

    /// <summary>Sample positions whose name is shared with another sample position.</summary>
    internal long NonUnique { get; set; }

    public void Dispose()
    {
        foreach (LeafResult l in Leaves) l.Dispose();
        Leaves.Clear();
    }
}

/// <summary>LIFO work list with termination detection (depth-first order keeps the number of pending tasks small).</summary>
internal sealed class WorkPool<T>
    where T : class
{
    private readonly Stack<T> _stack = new();
    private int _active;
    private Exception? _error;

    internal Exception? Error => _error;

    internal void Push(T t)
    {
        lock (_stack)
        {
            _stack.Push(t);
            Monitor.Pulse(_stack);
        }
    }

    /// <summary>Blocks until a task is available (true; the caller must call <see cref="Done"/>) or all work is finished or failed (false).</summary>
    internal bool Take(out T task)
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
/// The naming stage of PH1 (ADR-004 section 4.2): the sample positions are sorted by their X-prefix in <em>sample buckets</em> and equal X-prefixes get the same name. A bucket is a range
/// <c>[lo, hi)</c> of the total order (64-bit key, X-prefix, position) between two sample records; one task per bucket scans all sample positions in text order and keeps the records in its range.
/// A bucket of at most K_s records is sorted in RAM (radix on the key, the comparator on key ties), named, and written as a run sorted by slot. A bucket that holds more is split: one more scan counts its
/// records and keeps a uniform reservoir, whose quantiles are the splitters of the children (the root is such a bucket). Because the order includes the position, equal X-prefixes (a^n, satellites,
/// N runs) are spread over buckets like everything else, and no bucket exceeds K_s whatever the text. Equal groups that meet at a bucket boundary share one name (a boundary check after all buckets are known).
/// The names are the ranks of the groups in the order, starting at 1, as in <see cref="DcsBuilder"/>.
/// </summary>
internal sealed class SampleNamer
{
    private sealed record NamingTask(NamingBound? Lo, NamingBound? Hi, RangeNode Node, bool Sample, ulong Id);

    private readonly PackedText _text;
    private readonly SampleGeometry _g;
    private readonly ExtDcsOptions _o;
    private readonly ExtDcsLevelStats _stats;
    private readonly ArenaBudget _budget;
    private readonly string _tag;
    private readonly PrefixOrder _prefix;
    private readonly WorkPool<NamingTask> _pool = new();
    private readonly object _gate = new();
    private long _leaves;
    private long _scans;
    private long _splitScans;
    private long _fileRuns;

    internal SampleNamer(PackedText text, SampleGeometry geometry, ExtDcsOptions options, ExtDcsLevelStats stats, ArenaBudget budget, string tag)
    {
        _text = text;
        _g = geometry;
        _o = options;
        _stats = stats;
        _budget = budget;
        _tag = tag;
        _prefix = new PrefixOrder(text, geometry.Period);
    }

    /// <summary>Sorts and names the whole sample.</summary>
    internal NamingResult Run()
    {
        var root = new RangeNode();
        _pool.Push(new NamingTask(null, null, root, _g.Samples > _o.BucketCapacity, 1));
        int threads = Math.Max(1, _o.Threads);
        var extra = new List<Thread>();
        for (int t = 1; t < threads; t++)
        {
            var th = new Thread(Worker) { IsBackground = true, Name = $"ph1-naming-{t}" };
            extra.Add(th);
            th.Start();
        }
        Worker();
        foreach (Thread th in extra) th.Join();
        if (_pool.Error is { } error) throw new InvalidOperationException("sample naming failed: " + error.Message, error);

        var result = new NamingResult();
        var stack = new Stack<RangeNode>();
        stack.Push(root);
        while (stack.Count > 0)
        {
            RangeNode node = stack.Pop();
            if (node.Children is { } kids)
            {
                for (int i = kids.Count - 1; i >= 0; i--) stack.Push(kids[i]);
            }
            else if (node.Leaf is { Count: > 0 } leaf)
            {
                result.Leaves.Add(leaf);
            }
        }
        AssignNames(result);
        _stats.Leaves = _leaves;
        _stats.Scans = _scans;
        _stats.SplitScans = _splitScans;
        _stats.FileRuns = _fileRuns;
        _stats.DistinctNames = result.Distinct;
        _stats.NonUniqueRecords = result.NonUnique;
        return result;
    }

    /// <summary>
    /// Global names: walk the buckets in order, a group that continues over a bucket boundary keeps its name; also counts the groups of size one (the sample positions with a unique name).
    /// </summary>
    private void AssignNames(NamingResult r)
    {
        long names = 0;
        long openSize = 0;
        long singles = 0;
        LeafResult? prev = null;
        foreach (LeafResult leaf in r.Leaves)
        {
            bool continues = prev is not null && _text.Key(prev.LastPos) == _text.Key(leaf.FirstPos) && _prefix.Compare(prev.LastPos, leaf.FirstPos) == 0;
            long inner = 0;
            if (leaf.Distinct >= 2)
            {
                inner = leaf.Singles;
                if (leaf.FirstLen == 1) inner--;
                if (leaf.LastLen == 1) inner--;
            }
            if (continues)
            {
                leaf.Offset = (uint)(names - 1);
                names += leaf.Distinct - 1;
                if (leaf.Distinct == 1)
                {
                    openSize += leaf.Count;
                }
                else
                {
                    singles += inner;
                    openSize = leaf.LastLen;
                }
            }
            else
            {
                leaf.Offset = (uint)names;
                names += leaf.Distinct;
                if (openSize == 1) singles++;
                if (leaf.Distinct == 1)
                {
                    openSize = leaf.Count;
                }
                else
                {
                    if (leaf.FirstLen == 1) singles++;
                    singles += inner;
                    openSize = leaf.LastLen;
                }
            }
            prev = leaf;
        }
        if (openSize == 1) singles++;
        r.Distinct = names;
        r.NonUnique = _g.Samples - singles;
    }

    private void Worker()
    {
        using var arena = new BudgetedArena(_budget, null, _o.TouchBuffers);
        var w = new NamingWorker(this, arena);
        while (_pool.Take(out NamingTask? task))
        {
            try
            {
                w.Execute(task);
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
        lock (_gate)
        {
            _leaves += w.Leaves;
            _scans += w.Scans;
            _splitScans += w.SplitScans;
            _fileRuns += w.FileRuns;
        }
    }

    /// <summary>The merge of the runs into the name sink, one window of <c>ChunkSlots</c> slots at a time (a scatter by slot: no heap, no comparison).</summary>
    internal void Merge(NamingResult r, INameSink sink)
    {
        long s = _g.Samples;
        if (s == 0) return;
        using var arena = new BudgetedArena(_budget, null, _o.TouchBuffers);
        int windowSlots = (int)Math.Min(Math.Max(1, _o.ChunkSlots), s);
        Buf<uint> window = arena.Allocate<uint>(windowSlots, zero: false, "merge.window");
        int fileRuns = r.Leaves.Count(l => !l.Run!.InMemory);
        long fair = _o.MergeBufferBudget / 8 / Math.Max(1, fileRuns);
        int bufRecords = (int)Math.Clamp(fair, 1, Math.Max(1, _o.IoBufferBytes / 8));
        var readers = new List<(RunReader Reader, uint Offset)>(r.Leaves.Count);
        foreach (LeafResult leaf in r.Leaves)
        {
            Buf<ulong> buf = leaf.Run!.InMemory ? default : arena.Allocate<ulong>(bufRecords, zero: false, "merge.run");
            readers.Add((new RunReader(leaf.Run, buf, bufRecords), leaf.Offset));
        }
        for (long start = 0; start < s; start += windowSlots)
        {
            long end = Math.Min(s, start + windowSlots);
            Span<uint> span = window.AsSpan((int)(end - start));
            long got = 0;
            foreach ((RunReader reader, uint offset) in readers) got += reader.Drain((uint)end, offset, span, start);
            if (got != end - start) throw new InvalidOperationException($"merge window [{start}, {end}) received {got} records");
            sink.Accept(start, span);
        }
        long scattered = 0;
        foreach ((RunReader reader, _) in readers) scattered += reader.Scattered;
        if (scattered != s) throw new InvalidOperationException($"the runs hold {scattered} records for {s} sample positions");
    }

    /// <summary>The order of two sample records (key, X-prefix, position); for the splitter sample.</summary>
    private readonly struct RecordComparer(PackedText text, PrefixOrder prefix) : IComparer<uint>
    {
        public int Compare(uint a, uint b)
        {
            ulong ka = text.Key(a), kb = text.Key(b);
            if (ka != kb) return ka < kb ? -1 : 1;
            int c = prefix.Compare(a, b);
            return c != 0 ? c : a.CompareTo(b);
        }
    }

    /// <summary>The state of one worker thread: its arena, its bucket buffers (allocated at the first bucket), its counters.</summary>
    private sealed class NamingWorker(SampleNamer owner, BudgetedArena arena)
    {
        private Buf<ulong> _keys;
        private Buf<ulong> _keys2;
        private Buf<uint> _pos;
        private Buf<uint> _pos2;
        private bool _hasLeafBuffers;
        private uint[]? _reservoir;

        internal long Leaves;
        internal long Scans;
        internal long SplitScans;
        internal long FileRuns;

        internal void Execute(NamingTask task)
        {
            Scans++;
            if (task.Sample) Split(task);
            else Collect(task);
        }

        /// <summary>True when the record (key, p) lies in <c>[lo, hi)</c> of the order (key, X-prefix, position).</summary>
        private bool InRange(NamingBound? lo, NamingBound? hi, ulong key, uint p)
        {
            if (lo is not null)
            {
                if (key < lo.Key) return false;
                if (key == lo.Key && Cmp(lo.Pos, p) > 0) return false;
            }
            if (hi is not null)
            {
                if (key > hi.Key) return false;
                if (key == hi.Key && Cmp(p, hi.Pos) >= 0) return false;
            }
            return true;
        }

        /// <summary>Order of two records with equal keys: X-prefix, then position.</summary>
        private int Cmp(uint a, uint b)
        {
            int c = owner._prefix.Compare(a, b);
            return c != 0 ? c : a.CompareTo(b);
        }

        /// <summary>One scan that stores the records of the range in the bucket buffers; a bucket that outgrows K_s is queued again as a split.</summary>
        private void Collect(NamingTask task)
        {
            int cap = (int)Math.Min(owner._o.BucketCapacity, Math.Max(1, owner._g.Samples)); // a text with fewer sample positions than K_s needs no larger buffers
            EnsureLeafBuffers(cap);
            Span<ulong> keys = _keys.AsSpan(cap);
            Span<uint> pos = _pos.AsSpan(cap);
            SampleGeometry g = owner._g;
            PackedText text = owner._text;
            long n = g.N;
            int x = g.Period;
            int m = g.Classes;
            int[] residues = g.Cover.Residues;
            int stored = 0;
            bool overflow = false;
            for (long block = 0, blocks = n / x; block <= blocks && !overflow; block++)
            {
                long baseP = block * x;
                for (int c = 0; c < m; c++)
                {
                    long p = baseP + residues[c];
                    if (p > n) break;
                    ulong key = text.Key(p);
                    if (!InRange(task.Lo, task.Hi, key, (uint)p)) continue;
                    if (stored == cap)
                    {
                        overflow = true;
                        break;
                    }
                    keys[stored] = key;
                    pos[stored] = (uint)p;
                    stored++;
                }
            }
            if (overflow)
            {
                owner._pool.Push(task with { Sample = true });
                return;
            }
            Leaves++;
            task.Node.Leaf = Finish(stored);
        }

        /// <summary>One scan that counts the records of the range and keeps a uniform reservoir of them, then splits the range at the quantiles of the reservoir.</summary>
        private void Split(NamingTask task)
        {
            SplitScans++;
            int cap = Math.Max(3, owner._o.ReservoirSize);
            uint[] reservoir = _reservoir ??= new uint[cap];
            var rng = new Harness.SplitMix64(owner._o.Seed ^ (task.Id * 0x9E3779B97F4A7C15UL));
            SampleGeometry g = owner._g;
            PackedText text = owner._text;
            long n = g.N;
            int x = g.Period;
            int m = g.Classes;
            int[] residues = g.Cover.Residues;
            long count = 0;
            for (long block = 0, blocks = n / x; block <= blocks; block++)
            {
                long baseP = block * x;
                for (int c = 0; c < m; c++)
                {
                    long p = baseP + residues[c];
                    if (p > n) break;
                    ulong key = text.Key(p);
                    if (!InRange(task.Lo, task.Hi, key, (uint)p)) continue;
                    count++;
                    if (count <= cap)
                    {
                        reservoir[count - 1] = (uint)p;
                    }
                    else
                    {
                        ulong j = rng.Next() % (ulong)count;
                        if (j < (ulong)cap) reservoir[(int)j] = (uint)p;
                    }
                }
            }
            int size = (int)Math.Min(count, cap);
            if (size < 3) throw new InvalidOperationException($"a bucket of {count} records was queued for a split");
            Array.Sort(reservoir, 0, size, new RecordComparer(text, owner._prefix));
            long bucket = owner._o.BucketCapacity;
            int fan = (int)Math.Min(Math.Max(2, (2 * count + bucket - 1) / bucket), Math.Max(2, owner._o.FanOut));
            var kids = new List<RangeNode>(fan);
            var tasks = new List<NamingTask>(fan);
            NamingBound? prev = task.Lo;
            int last = 0;
            ulong id = task.Id * 0x100000001B3UL;
            for (int t = 1; t < fan; t++)
            {
                int i = (int)((long)t * size / fan);
                if (i <= last) continue; // at least one reservoir record lies before the splitter, and the splitters strictly increase
                last = i;
                var bound = new NamingBound(reservoir[i], text.Key(reservoir[i]));
                var node = new RangeNode();
                kids.Add(node);
                tasks.Add(new NamingTask(prev, bound, node, false, id + (ulong)kids.Count));
                prev = bound;
            }
            var lastNode = new RangeNode();
            kids.Add(lastNode);
            tasks.Add(new NamingTask(prev, task.Hi, lastNode, false, id + (ulong)kids.Count));
            task.Node.Children = kids;
            for (int i = tasks.Count - 1; i >= 0; i--) owner._pool.Push(tasks[i]);
        }

        private void EnsureLeafBuffers(int cap)
        {
            if (_hasLeafBuffers) return;
            _keys = arena.Allocate<ulong>(cap, zero: false, "bucket.keys");
            _keys2 = arena.Allocate<ulong>(cap, zero: false, "bucket.keys2");
            _pos = arena.Allocate<uint>(cap, zero: false, "bucket.pos");
            _pos2 = arena.Allocate<uint>(cap, zero: false, "bucket.pos2");
            _hasLeafBuffers = true;
        }

        /// <summary>Sorts the <paramref name="m"/> collected records by (key, X-prefix, position), names the groups, and writes the run sorted by slot.</summary>
        private LeafResult Finish(int m)
        {
            if (m == 0) return new LeafResult();
            Span<ulong> keys = _keys.AsSpan(m);
            Span<uint> pos = _pos.AsSpan(m);
            if (m <= 24) InsertionSort(keys, pos);
            else RadixSort(m);
            keys = _keys.AsSpan(m);
            pos = _pos.AsSpan(m);
            int s = 0;
            while (s < m)
            {
                int e = s + 1;
                while (e < m && keys[e] == keys[s]) e++;
                if (e - s >= 2) MergeSort(pos.Slice(s, e - s), _pos2.AsSpan(e - s));
                s = e;
            }

            Span<ulong> values = _keys2.AsSpan(m);
            long distinct = 0, firstLen = 0, singles = 0;
            int groupStart = 0;
            for (int i = 0; i < m; i++)
            {
                if (i == 0 || keys[i] != keys[i - 1] || owner._prefix.Compare(pos[i - 1], pos[i]) != 0)
                {
                    if (i > 0)
                    {
                        int size = i - groupStart;
                        if (size == 1) singles++;
                        if (distinct == 1) firstLen = size;
                    }
                    distinct++;
                    groupStart = i;
                }
                values[i] = ((ulong)owner._g.SlotOf(pos[i]) << 32) | (uint)distinct;
            }
            long lastLen = m - groupStart;
            if (lastLen == 1) singles++;
            if (distinct == 1) firstLen = lastLen;
            uint firstPos = pos[0], lastPos = pos[m - 1];
            values.Sort();
            BucketRun run = BucketRun.Create(values, owner._o.ScratchDir, owner._tag, owner._o.RunSpillBytes);
            if (!run.InMemory) FileRuns++;
            return new LeafResult { Count = m, Distinct = distinct, FirstPos = firstPos, LastPos = lastPos, FirstLen = firstLen, LastLen = lastLen, Singles = singles, Run = run };
        }

        /// <summary>
        /// Stable merge sort of a run of positions with equal keys by X-prefix (binary insertion up to six). The positions arrive in increasing order, so equal prefixes keep that order, which is the
        /// position tie-break of the bucket order; a merge is skipped when the halves are already in order (a^n, tandem repeats).
        /// </summary>
        private void MergeSort(Span<uint> a, Span<uint> tmp)
        {
            int n = a.Length;
            PrefixOrder prefix = owner._prefix;
            if (n <= 6)
            {
                for (int i = 1; i < n; i++)
                {
                    uint v = a[i];
                    int lo = 0, hi = i;
                    while (lo < hi)
                    {
                        int mid = (lo + hi) >> 1;
                        if (prefix.Compare(a[mid], v) <= 0) lo = mid + 1;
                        else hi = mid;
                    }
                    for (int j = i; j > lo; j--) a[j] = a[j - 1];
                    a[lo] = v;
                }
                return;
            }
            int half = n >> 1;
            MergeSort(a[..half], tmp[..half]);
            MergeSort(a[half..], tmp[half..]);
            if (prefix.Compare(a[half - 1], a[half]) <= 0) return;
            a[..half].CopyTo(tmp);
            int p = 0, q = half, k = 0;
            while (p < half && q < n)
            {
                if (prefix.Compare(a[q], tmp[p]) < 0) a[k++] = a[q++];
                else a[k++] = tmp[p++];
            }
            while (p < half) a[k++] = tmp[p++];
        }

        private static void InsertionSort(Span<ulong> keys, Span<uint> pos)
        {
            for (int i = 1; i < keys.Length; i++)
            {
                ulong k = keys[i];
                uint p = pos[i];
                int j = i - 1;
                while (j >= 0 && keys[j] > k)
                {
                    keys[j + 1] = keys[j];
                    pos[j + 1] = pos[j];
                    j--;
                }
                keys[j + 1] = k;
                pos[j + 1] = p;
            }
        }

        /// <summary>LSD radix sort of (key, position) on the 64-bit key, stable, skipping the digits on which all keys agree.</summary>
        private void RadixSort(int m)
        {
            Span<ulong> k1 = _keys.AsSpan(m);
            Span<ulong> k2 = _keys2.AsSpan(m);
            Span<uint> p1 = _pos.AsSpan(m);
            Span<uint> p2 = _pos2.AsSpan(m);
            Span<int> hist = stackalloc int[8 * 256];
            hist.Clear();
            foreach (ulong k in k1)
            {
                hist[(int)(k & 0xFF)]++;
                hist[256 + (int)((k >> 8) & 0xFF)]++;
                hist[512 + (int)((k >> 16) & 0xFF)]++;
                hist[768 + (int)((k >> 24) & 0xFF)]++;
                hist[1024 + (int)((k >> 32) & 0xFF)]++;
                hist[1280 + (int)((k >> 40) & 0xFF)]++;
                hist[1536 + (int)((k >> 48) & 0xFF)]++;
                hist[1792 + (int)(k >> 56)]++;
            }
            ulong first = k1[0];
            bool swapped = false;
            for (int d = 0; d < 8; d++)
            {
                Span<int> h = hist.Slice(d * 256, 256);
                int shift = d * 8;
                if (h[(int)((first >> shift) & 0xFF)] == m) continue;
                int sum = 0;
                for (int b = 0; b < 256; b++)
                {
                    int c = h[b];
                    h[b] = sum;
                    sum += c;
                }
                Span<ulong> sk = swapped ? k2 : k1;
                Span<ulong> dk = swapped ? k1 : k2;
                Span<uint> sp = swapped ? p2 : p1;
                Span<uint> dp = swapped ? p1 : p2;
                for (int i = 0; i < m; i++)
                {
                    ulong k = sk[i];
                    int at = h[(int)((k >> shift) & 0xFF)]++;
                    dk[at] = k;
                    dp[at] = sp[i];
                }
                swapped = !swapped;
            }
            if (swapped)
            {
                k2.CopyTo(k1);
                p2.CopyTo(p1);
            }
        }
    }
}
