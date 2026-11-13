using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph4;

namespace SuffixTree.EmProto.Ph5;

/// <summary>How the single-array child table (layout A) encodes a field: <see cref="Tagged"/> = kind in bits 30-31 (the in-RAM <c>TaggedChild</c>), <see cref="Packed"/> = the raw value (the in-RAM <c>PackedChild</c>).</summary>
internal enum ChildCodec
{
    Tagged,
    Packed,
}

/// <summary>Parameters of PH5 (ADR-004 section 3.2). The output never depends on any of them.</summary>
internal sealed record ChildBuildOptions
{
    internal required string ScratchDir { get; init; }

    /// <summary>Entries per stack block (ADR: 256 K); the two top blocks are in RAM.</summary>
    internal int BlockEntries { get; init; } = 256 * 1024;

    /// <summary>Size of one mapped window of the CHILD, PSV and NSV writers (a power of two).</summary>
    internal int WindowBytes { get; init; } = 1 << 20;

    /// <summary>Windows the layout A writer keeps mapped at a time.</summary>
    internal int WindowSlots { get; init; } = 8;

    /// <summary>msync a window when it is released (bounds the dirty pages).</summary>
    internal bool FlushWindows { get; init; } = true;

    /// <summary>fsync the spill file after this many unflushed bytes (0 = never).</summary>
    internal long SpillFlushBytes { get; init; } = 16L << 20;

    internal int ReadBufferBytes { get; init; } = 1 << 20;

    internal ChildCodec Codec { get; init; } = ChildCodec.Tagged;

    internal ArenaBudget? Budget { get; init; }

    internal ArenaOptions? ArenaGuards { get; init; }

    internal bool TouchBuffers { get; init; }

    /// <summary>Layout F from a file: run the backward NSV scan on its own thread while the forward scan runs (the two passes are independent: both read the LCP file, they write different files); false = one pass after the other.</summary>
    internal bool ParallelBackward { get; init; }
}

/// <summary>The files the forward scan writes (null = that part is off): the single-array child table (layout A) and the PSV array of layout F.</summary>
internal sealed record ChildScanOutputs
{
    internal string? ChildPath { get; init; }

    internal string? PsvPath { get; init; }
}

/// <summary>What the forward scan found and did.</summary>
internal sealed class ChildScanResult
{
    internal long N { get; set; }

    /// <summary>Internal nodes with a value above 0 (the root is not counted).</summary>
    internal long InternalCount { get; set; }

    /// <summary>NODE_COUNT = <c>1 + I + n + 1</c> (root, internal nodes, one leaf per suffix including the empty one); 1 for the empty text.</summary>
    internal long NodeCount { get; set; }

    internal uint LrsDepth { get; set; }

    /// <summary>Left boundary of the first-closing deepest interval: the longest repeated substring starts at <c>SA[LrsRank]</c>.</summary>
    internal uint LrsRank { get; set; }

    internal long MaxLcp { get; set; }

    /// <summary>Peak height of the scan stack (at most <c>1 + max LCP</c>).</summary>
    internal long StackPeak { get; set; }

    internal SpillStats Spill { get; set; }

    internal WindowStats ChildWindows { get; set; }

    internal WindowStats PsvWindows { get; set; }

    /// <summary>Seconds spent inside <see cref="ChildScan.Consume"/> and <see cref="ChildScan.Complete"/> (the scan itself, without reading the LCP).</summary>
    internal double Seconds { get; set; }
}

/// <summary>An open interval of the scan: value, left boundary, first and last l-index (<see cref="ChildScan.None"/> = none yet). 16 bytes (ADR-004 section 3.2: "12-16 B").</summary>
[StructLayout(LayoutKind.Sequential)]
internal struct IntervalEntry
{
    internal uint Ell;
    internal uint Lb;
    internal uint First;
    internal uint Last;
}

/// <summary>Events of the forward scan that the writers of layout A and layout F listen to (the in-RAM <c>IIntervalSink</c> plus the position event of the PSV).</summary>
internal interface IChildSink
{
    /// <summary><paramref name="k"/> is an l-index of the open interval; <paramref name="previous"/> is its previous l-index (the link <c>next[previous] = k</c>) or <see cref="ChildScan.None"/> for the first.</summary>
    void LIndex(uint previous, uint k);

    /// <summary>The internal interval <c>[lb..rb]</c> with value <paramref name="ell"/> and first l-index <paramref name="first"/> closes (children before their parent, the root last); <paramref name="lastChild"/>: its right boundary is its parent's.</summary>
    void Close(uint ell, uint lb, uint rb, uint first, bool lastChild);

    /// <summary>Rank <paramref name="i"/> has been processed; <paramref name="psv"/> is its previous smaller value (strict) or <see cref="ChildScan.None"/>.</summary>
    void Position(uint i, uint psv);
}

internal struct NoSink : IChildSink
{
    public readonly void LIndex(uint previous, uint k)
    {
    }

    public readonly void Close(uint ell, uint lb, uint rb, uint first, bool lastChild)
    {
    }

    public readonly void Position(uint i, uint psv)
    {
    }
}

internal struct TeeSink<TA, TB>(TA a, TB b) : IChildSink
    where TA : struct, IChildSink
    where TB : struct, IChildSink
{
    private TA _a = a;
    private TB _b = b;

    public void LIndex(uint previous, uint k)
    {
        _a.LIndex(previous, k);
        _b.LIndex(previous, k);
    }

    public void Close(uint ell, uint lb, uint rb, uint first, bool lastChild)
    {
        _a.Close(ell, lb, rb, first, lastChild);
        _b.Close(ell, lb, rb, first, lastChild);
    }

    public void Position(uint i, uint psv)
    {
        _a.Position(i, psv);
        _b.Position(i, psv);
    }
}

/// <summary>Layout A, tagged codec: the writes of the in-RAM <c>TaggedChild</c> (next link at the previous l-index, up at <c>rb</c>, down at <c>lb</c>), kind in bits 30-31.</summary>
internal readonly struct TaggedChildSink(AlignedWindowWriter writer) : IChildSink
{
    private const uint NextTag = 1u << 30;
    private const uint UpTag = 2u << 30;
    private const uint DownTag = 3u << 30;

    public void LIndex(uint previous, uint k)
    {
        if (previous != ChildScan.None) writer.Write32(previous, NextTag | k);
    }

    public void Close(uint ell, uint lb, uint rb, uint first, bool lastChild)
    {
        if (lastChild) writer.Write32(lb, DownTag | first);
        else writer.Write32(rb, UpTag | first);
    }

    public readonly void Position(uint i, uint psv)
    {
    }
}

/// <summary>Layout A, packed codec: the writes of the in-RAM <c>PackedChild</c> (the raw value; a next link only when it is not <c>k + 1</c>).</summary>
internal readonly struct PackedChildSink(AlignedWindowWriter writer) : IChildSink
{
    public void LIndex(uint previous, uint k)
    {
        if (previous != ChildScan.None && k != previous + 1) writer.Write32(previous, k);
    }

    public void Close(uint ell, uint lb, uint rb, uint first, bool lastChild) => writer.Write32(lastChild ? lb : rb, first);

    public readonly void Position(uint i, uint psv)
    {
    }
}

/// <summary>Layout F, forward half: PSV[i] written at entry i, strictly ascending (the writes are sequential).</summary>
internal readonly struct PsvSink(AlignedWindowWriter writer) : IChildSink
{
    public readonly void LIndex(uint previous, uint k)
    {
    }

    public readonly void Close(uint ell, uint lb, uint rb, uint first, bool lastChild)
    {
    }

    public void Position(uint i, uint psv) => writer.Write32(i, psv);
}

/// <summary>The scan behind <see cref="ChildScan"/>, specialised by the sink (the JIT inlines the writers).</summary>
internal interface IScanCore : IDisposable
{
    void Consume(long firstRank, ReadOnlySpan<uint> lcp);

    ChildScanResult Complete();
}

internal sealed class ScanCore<TSink> : IScanCore
    where TSink : struct, IChildSink
{
    private readonly long _n;
    private readonly BudgetedArena _arena;
    private readonly SpillStack<IntervalEntry> _stack;
    private readonly AlignedWindowWriter? _child;
    private readonly AlignedWindowWriter? _psv;
    private TSink _sink;
    private long _nextRank;
    private long _internal;
    private uint _lrsDepth;
    private uint _lrsRank;
    private uint _maxLcp;
    private long _ticks;
    private bool _completed;
    private bool _disposed;

    internal ScanCore(long n, TSink sink, AlignedWindowWriter? child, AlignedWindowWriter? psv, BudgetedArena arena, string spillPath, ChildBuildOptions options)
    {
        _n = n;
        _sink = sink;
        _child = child;
        _psv = psv;
        _arena = arena;
        _stack = new SpillStack<IntervalEntry>(arena, spillPath, options.BlockEntries, options.SpillFlushBytes);
        _stack.Push(new IntervalEntry { Ell = 0, Lb = 0, First = ChildScan.None, Last = ChildScan.None });
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public void Consume(long firstRank, ReadOnlySpan<uint> lcp)
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        if (firstRank != _nextRank) throw new InvalidOperationException($"LCP values must arrive in rank order: expected rank {_nextRank}, got {firstRank}");
        if (firstRank + lcp.Length > _n + 1) throw new InvalidOperationException($"ranks {firstRank}..{firstRank + lcp.Length - 1} exceed the {_n + 1} ranks of a text of {_n} symbols");
        long t0 = Stopwatch.GetTimestamp();
        for (int j = 0; j < lcp.Length; j++)
        {
            long rank = firstRank + j;
            uint cur = lcp[j];
            if (rank <= 1)
            {
                if (cur != 0) throw new InvalidDataException($"LCP[{rank}] = {cur}: the values at ranks 0 and 1 must be 0");
                if (rank == 0)
                {
                    _sink.Position(0, ChildScan.None);
                    continue;
                }
            }
            else if (cur > _n)
            {
                throw new InvalidDataException($"LCP[{rank}] = {cur} exceeds the text length {_n}");
            }
            if (cur > _maxLcp) _maxLcp = cur;
            Step((uint)rank, cur);
        }
        _nextRank = firstRank + lcp.Length;
        _ticks += Stopwatch.GetTimestamp() - t0;
    }

    public ChildScanResult Complete()
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        if (_completed) throw new InvalidOperationException("the scan is already complete");
        if (_nextRank != _n + 1) throw new InvalidOperationException($"{_nextRank} of the {_n + 1} LCP values were delivered");
        long t0 = Stopwatch.GetTimestamp();
        if (_n > 0)
        {
            while (_stack.Height > 0)
            {
                IntervalEntry e = _stack.Top;
                _stack.Pop();
                Close(e.Ell, e.Lb, (uint)_n, e.First, true);
            }
        }
        _completed = true;
        _child?.Dispose();
        _psv?.Dispose();
        _ticks += Stopwatch.GetTimestamp() - t0;
        return new ChildScanResult
        {
            N = _n,
            InternalCount = _internal,
            NodeCount = _n == 0 ? 1 : 1 + _internal + _n + 1,
            LrsDepth = _lrsDepth,
            LrsRank = _lrsRank,
            MaxLcp = _maxLcp,
            StackPeak = _stack.PeakHeight,
            Spill = _stack.Stats,
            ChildWindows = _child?.Stats ?? default,
            PsvWindows = _psv?.Stats ?? default,
            Seconds = _ticks / (double)Stopwatch.Frequency,
        };
    }

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        _child?.Dispose();
        _psv?.Dispose();
        _stack.Dispose();
        _arena.Dispose();
    }

    /// <summary>One rank of the bottom-up scan of the lcp-interval tree (the algorithm of the in-RAM <c>LcpIntervalScan</c>, with the stack in a <see cref="SpillStack{T}"/>).</summary>
    private void Step(uint i, uint cur)
    {
        uint left = i - 1;
        ref IntervalEntry top = ref _stack.Top;
        while (top.Ell > cur)
        {
            IntervalEntry e = top;
            _stack.Pop();
            top = ref _stack.Top;
            Close(e.Ell, e.Lb, i - 1, e.First, cur < top.Ell);
            left = e.Lb;
        }
        uint psv;
        if (top.Ell < cur)
        {
            _stack.Push(new IntervalEntry { Ell = cur, Lb = left, First = i, Last = i });
            _sink.LIndex(ChildScan.None, i);
            psv = left;
        }
        else if (top.First == ChildScan.None)
        {
            ref IntervalEntry root = ref _stack.TopForUpdate();
            root.First = i;
            root.Last = i;
            _sink.LIndex(ChildScan.None, i);
            psv = ChildScan.None;
        }
        else
        {
            psv = top.Ell == 0 ? ChildScan.None : top.Lb;
            ref IntervalEntry t = ref _stack.TopForUpdate();
            _sink.LIndex(t.Last, i);
            t.Last = i;
        }
        _sink.Position(i, psv);
    }

    private void Close(uint ell, uint lb, uint rb, uint first, bool lastChild)
    {
        if (ell > 0)
        {
            _internal++;
            if (ell > _lrsDepth)
            {
                _lrsDepth = ell;
                _lrsRank = lb;
            }
        }
        _sink.Close(ell, lb, rb, first, lastChild);
    }
}

/// <summary>
/// The forward scan of PH5 (ADR-004 sections 3.2 and 4.2, row PH5): one pass over the LCP in rank order with the external stack, writing layout A (the child table in the single-array form, writes behind the cursor through
/// bounded mapped windows) and the PSV half of layout F (sequential), and counting the internal nodes and the deepest interval (NODE_COUNT, LRS_DEPTH, LRS_RANK). It is an <see cref="ILcpChunkConsumer"/>: the PH4 hook
/// (<c>Consume(firstRank, values)</c>, dense resolved values in ascending rank order) feeds it fused with PH4-B3, and <see cref="ExternalChild"/> feeds it from a compact LCP file.
/// </summary>
internal sealed class ChildScan : ILcpChunkConsumer, IDisposable
{
    /// <summary>"No index": the stored value of a missing previous smaller value (-1) and of a missing l-index.</summary>
    internal const uint None = uint.MaxValue;

    private readonly IScanCore _core;

    private ChildScan(IScanCore core) => _core = core;

    /// <summary>Creates the scan of a text of <paramref name="n"/> symbols; the output files are created at their final length (<c>4 x (n + 1)</c> bytes) at once.</summary>
    /// <exception cref="NotSupportedException">Tagged codec for a text above 2^30 - 1 symbols, or a text of 2^32 - 1 symbols and more (32-bit positions).</exception>
    internal static ChildScan Create(long n, ChildScanOutputs outputs, ChildBuildOptions options)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        ArgumentNullException.ThrowIfNull(outputs);
        ArgumentNullException.ThrowIfNull(options);
        if (n >= uint.MaxValue - 1) throw new NotSupportedException("positions are 32-bit in the prototype");
        if (outputs.ChildPath is not null && options.Codec == ChildCodec.Tagged && n > TaggedChild.MaxLength) throw new NotSupportedException($"the tagged codec holds ranks in 30 bits; a text of {n} symbols needs the packed codec");
        Directory.CreateDirectory(options.ScratchDir);
        var arena = new BudgetedArena(options.Budget ?? new ArenaBudget(long.MaxValue), options.ArenaGuards, options.TouchBuffers);
        AlignedWindowWriter? child = null, psv = null;
        try
        {
            if (outputs.ChildPath is not null) child = new AlignedWindowWriter(outputs.ChildPath, 4 * (n + 1), options.WindowBytes, options.WindowSlots, options.FlushWindows);
            if (outputs.PsvPath is not null) psv = new AlignedWindowWriter(outputs.PsvPath, 4 * (n + 1), options.WindowBytes, 2, options.FlushWindows);
            string spill = Path.Combine(options.ScratchDir, "ph5-stack-" + Guid.NewGuid().ToString("N") + ".spill");
            IScanCore core = (child, psv, options.Codec) switch
            {
                (null, null, _) => new ScanCore<NoSink>(n, default, null, null, arena, spill, options),
                (not null, null, ChildCodec.Tagged) => new ScanCore<TaggedChildSink>(n, new TaggedChildSink(child), child, null, arena, spill, options),
                (not null, null, _) => new ScanCore<PackedChildSink>(n, new PackedChildSink(child), child, null, arena, spill, options),
                (null, not null, _) => new ScanCore<PsvSink>(n, new PsvSink(psv), null, psv, arena, spill, options),
                (not null, not null, ChildCodec.Tagged) => new ScanCore<TeeSink<TaggedChildSink, PsvSink>>(n, new TeeSink<TaggedChildSink, PsvSink>(new TaggedChildSink(child), new PsvSink(psv)), child, psv, arena, spill, options),
                (not null, not null, _) => new ScanCore<TeeSink<PackedChildSink, PsvSink>>(n, new TeeSink<PackedChildSink, PsvSink>(new PackedChildSink(child), new PsvSink(psv)), child, psv, arena, spill, options),
            };
            return new ChildScan(core);
        }
        catch
        {
            child?.Dispose();
            psv?.Dispose();
            arena.Dispose();
            throw;
        }
    }

    /// <summary>The next block of LCP values: ranks <c>firstRank ..</c>, which must continue where the previous call ended (the first call starts at rank 0).</summary>
    /// <exception cref="InvalidOperationException">The ranks are not in order or exceed n + 1.</exception>
    /// <exception cref="InvalidDataException">LCP[0] or LCP[1] is not 0, or a value exceeds the text length.</exception>
    public void Consume(long firstRank, ReadOnlySpan<uint> lcp) => _core.Consume(firstRank, lcp);

    /// <summary>Closes the intervals still open, flushes and closes the output files and returns what the scan found.</summary>
    /// <exception cref="InvalidOperationException">Not all <c>n + 1</c> values were delivered.</exception>
    internal ChildScanResult Complete() => _core.Complete();

    public void Dispose() => _core.Dispose();
}
