namespace SuffixTree.Esa;

/// <summary>
/// The write pattern of the child-table scan, replayed per layout from the scan events (prototype P3; ADR-004 section 3.2 requires the child-table writes
/// of a streaming builder to be sequential in one direction or LIFO-local). The scan cursor is the rank the scan has just read (<c>k</c> for an l-index,
/// <c>rb + 1</c> for a closing interval). A write is AT THE CURSOR when it targets rank <c>cursor - 1</c> or <c>cursor</c> (an append in the forward direction)
/// and BEHIND it otherwise (the left boundary of a closing interval, the previous l-index of an open interval). Page misses are those of an LRU cache of
/// <c>cachePages</c> pages of <c>pageSlots</c> four-byte slots (4,096-byte pages: 1,024 slots), counted per array; a miss is a page that a mapped file would have to
/// read in before the write (and write back later).
/// </summary>
internal sealed class ChildWriteTrace
{
    private readonly bool _trackSlots;
    private readonly int _pageSlots;
    private readonly int _cachePages;
    private readonly HashSet<long>? _slots;
    private readonly Dictionary<long, long> _pages = [];
    private long _tick;
    private long _previousBehindTarget = long.MaxValue;
    private bool _previousWasBehind;

    private ChildWriteTrace(bool trackSlots, int pageSlots, int cachePages)
    {
        _trackSlots = trackSlots;
        _pageSlots = pageSlots;
        _cachePages = cachePages;
        if (trackSlots) _slots = [];
    }

    /// <summary>Total number of writes.</summary>
    internal long Writes { get; private set; }

    /// <summary>Writes at the cursor (an append in the forward direction).</summary>
    internal long AtCursor { get; private set; }

    /// <summary>Writes behind the cursor.</summary>
    internal long Behind { get; private set; }

    /// <summary>The largest distance <c>cursor - target</c> of a write behind the cursor.</summary>
    internal long MaxBackDistance { get; private set; }

    /// <summary>Number of maximal runs of consecutive behind-writes whose targets strictly descend (a^n: one run).</summary>
    internal long BackwardRuns { get; private set; }

    /// <summary>Page loads of the LRU cache over all writes.</summary>
    internal long PageMisses { get; private set; }

    /// <summary>Writes to a slot (array, rank) that was written before; always 0 when the slot is tracked and the layout is sound.</summary>
    internal long DuplicateSlots { get; private set; }

    /// <summary>Peak height of the scan stack.</summary>
    internal int StackPeak { get; private set; }

    /// <summary>Replays the scan of <paramref name="lcp"/> as the writes of <paramref name="layout"/>.</summary>
    internal static ChildWriteTrace Run(ChildLayoutKind layout, int n, CompactLcp lcp, int pageSlots = 1024, int cachePages = 64, bool trackSlots = true)
    {
        ArgumentNullException.ThrowIfNull(lcp);
        ArgumentOutOfRangeException.ThrowIfLessThan(pageSlots, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(cachePages, 1);
        var trace = new ChildWriteTrace(trackSlots, pageSlots, cachePages);
        var stream = new CompactLcpStream(lcp);
        var sink = new TraceSink(trace, layout);
        trace.StackPeak = LcpIntervalScan.Run(n, ref stream, ref sink);
        return trace;
    }

    private void Write(int array, long index, long cursor)
    {
        Writes++;
        if (_trackSlots && !_slots!.Add(((long)array << 40) | index)) DuplicateSlots++;
        if (index >= cursor - 1)
        {
            AtCursor++;
            _previousWasBehind = false;
        }
        else
        {
            Behind++;
            MaxBackDistance = Math.Max(MaxBackDistance, cursor - index);
            if (!_previousWasBehind || index >= _previousBehindTarget) BackwardRuns++;
            _previousWasBehind = true;
            _previousBehindTarget = index;
        }

        long page = ((long)array << 40) | (index / _pageSlots);
        _tick++;
        if (_pages.ContainsKey(page))
        {
            _pages[page] = _tick;
            return;
        }
        PageMisses++;
        if (_pages.Count >= _cachePages)
        {
            long oldestPage = 0, oldestTick = long.MaxValue;
            foreach (KeyValuePair<long, long> entry in _pages)
                if (entry.Value < oldestTick)
                {
                    oldestTick = entry.Value;
                    oldestPage = entry.Key;
                }
            _pages.Remove(oldestPage);
        }
        _pages[page] = _tick;
    }

    private readonly struct TraceSink(ChildWriteTrace trace, ChildLayoutKind layout) : IIntervalSink
    {
        public void LIndex(int previous, int k)
        {
            if (previous < 0) return;
            switch (layout)
            {
                case ChildLayoutKind.Reference3:
                    trace.Write(0, previous, k);
                    break;
                case ChildLayoutKind.Packed:
                    if (k != previous + 1) trace.Write(0, previous, k);
                    break;
                default:
                    trace.Write(0, previous, k);
                    break;
            }
        }

        public void Close(int ell, int lb, int rb, int first, bool lastChild)
        {
            long cursor = rb + 1L;
            if (layout == ChildLayoutKind.Reference3)
            {
                if (lastChild) trace.Write(2, lb, cursor);
                else trace.Write(1, rb, cursor);
            }
            else
            {
                trace.Write(0, lastChild ? lb : rb, cursor);
            }
        }
    }
}
