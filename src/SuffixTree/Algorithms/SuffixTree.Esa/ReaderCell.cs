using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace SuffixTree.Esa;

/// <summary>
/// The thread-exclusive slot of the readers of the file-backed trees (work package 1b.9g; the protocol and its proof are in ADR-004 section 5a, "1b.9g"). One cell belongs to one thread at a time. While the thread is inside a call of a file-backed
/// tree, the cell holds the id of that tree in <c>_first</c> (the outermost call) or in a segment of the nested calls; it holds 0 otherwise. Only the owner thread writes the entries, with plain stores (<see cref="Volatile.Write{T}(ref T, T)"/> is a plain
/// store on x64 and a store-release on arm64): no read-modify-write and no fence on the path of a call. The releaser of a disposed tree (<see cref="ReaderGate{TSeams}"/>) reads the entries of every cell after a process-wide barrier
/// (<see cref="Interlocked.MemoryBarrierProcessWide"/>), which is what makes the plain stores enough. The cells of all threads are in one list that is only appended to; a cell whose thread has ended is taken over by the next thread that registers, so the list is
/// as long as the largest number of threads that were alive at one time, whatever the number of trees and of short-lived threads. A cell is 320 bytes with the written fields in the middle (128 bytes of room on each side): two cells never share a cache line, and neither
/// does a cell with the neighbour that the allocator places after it.
/// </summary>
[StructLayout(LayoutKind.Explicit, Size = CellSize)]
internal sealed class ReaderCell
{
    private const int CellSize = 320;
    private const int HotOffset = 128;
    private const int SegmentLength = 8;

    [ThreadStatic]
    private static ReaderCell? t_current;

    private static readonly Lock s_lock = new();
    private static ReaderCell? s_head;
    private static int s_count;
    private static long s_lastId;

    // cold: written once, before the cell is in the list (the list link), or under the lock (the owner)
    [FieldOffset(0)]
    private ReaderCell? _next;

    [FieldOffset(8)]
    private Thread? _owner;

    // hot: written by the owner thread only
    [FieldOffset(HotOffset)]
    private long _first;

    [FieldOffset(HotOffset + 8)]
    private int _nestedTop;

    [FieldOffset(HotOffset + 16)]
    private Segment? _segments;

    private sealed class Segment
    {
        internal readonly long[] Ids = new long[SegmentLength];
        internal Segment? Next;
    }

    /// <summary>The cell of the calling thread (registered on the first use of the thread).</summary>
    internal static ReaderCell Current
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => t_current ?? Register();
    }

    /// <summary>A new id for a gate: never 0 (an empty cell holds 0) and never used twice in the life of the process.</summary>
    internal static long NextId() => Interlocked.Increment(ref s_lastId);

    /// <summary>The number of cells that exist (the largest number of threads that were registered at one time, not the number of threads that ever were).</summary>
    internal static int CellCount => Volatile.Read(ref s_count);

    /// <summary>The number of calls this cell is inside of right now (tests).</summary>
    internal int Depth => (_first == 0 ? 0 : 1) + _nestedTop;

    /// <summary>The calling thread enters a call of the tree <paramref name="id"/>: the entry is stored, nothing else (the caller reads the flag of the tree afterwards).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal void Publish(long id)
    {
        if (_first == 0) Volatile.Write(ref _first, id);
        else PublishNested(id);
    }

    /// <summary>The calling thread leaves its innermost call: the entry is cleared with a release store, so that everything the call read is read before the entry is gone.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal void Retract()
    {
        if (_nestedTop == 0) Volatile.Write(ref _first, 0);
        else RetractNested();
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    private void PublishNested(long id)
    {
        int index = _nestedTop;
        Segment segment = SegmentOf(index, create: true)!;
        Volatile.Write(ref segment.Ids[index % SegmentLength], id);
        _nestedTop = index + 1;
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    private void RetractNested()
    {
        int index = _nestedTop - 1;
        Segment segment = SegmentOf(index, create: false)!;
        Volatile.Write(ref segment.Ids[index % SegmentLength], 0);
        _nestedTop = index;
    }

    // the segment that holds the nested entry <index>; a new segment is linked (a release store) before an entry is stored in it, so that a scan that finds the entry has found the link
    private Segment? SegmentOf(int index, bool create)
    {
        ref Segment? link = ref _segments;
        for (int hops = index / SegmentLength; hops > 0; hops--) link = ref link!.Next;
        Segment? segment = link;
        if (segment is null && create)
        {
            segment = new Segment();
            Volatile.Write(ref link, segment);
        }

        return segment;
    }

    /// <summary>
    /// Whether any thread is inside a call of the tree <paramref name="id"/>: every cell of the list, read with acquire loads. Only meaningful after a process-wide barrier that follows the setting of the tree's flag
    /// (see <see cref="ReaderGate{TSeams}"/>): without it a thread whose entry is still in its store buffer would not be seen.
    /// </summary>
    internal static bool AnyInside(long id)
    {
        for (ReaderCell? cell = Volatile.Read(ref s_head); cell is not null; cell = cell._next)
            if (cell.Holds(id)) return true;
        return false;
    }

    private bool Holds(long id)
    {
        if (Volatile.Read(ref _first) == id) return true;
        for (Segment? segment = Volatile.Read(ref _segments); segment is not null; segment = Volatile.Read(ref segment.Next))
            for (int i = 0; i < SegmentLength; i++)
                if (Volatile.Read(ref segment.Ids[i]) == id) return true;
        return false;
    }

    /// <summary>The first call of a thread: a cell of a thread that has ended, or a new one, in the list before the thread's first entry.</summary>
    [MethodImpl(MethodImplOptions.NoInlining)]
    private static ReaderCell Register()
    {
        Thread me = Thread.CurrentThread;
        ReaderCell cell;
        lock (s_lock)
        {
            ReaderCell? ended = null;
            for (ReaderCell? c = s_head; c is not null && ended is null; c = c._next)
                if (c._owner is { IsAlive: false }) ended = c;
            if (ended is not null)
            {
                ended.Clear();
                ended._owner = me;
                cell = ended;
            }
            else
            {
                cell = new ReaderCell { _owner = me, _next = s_head };
                Volatile.Write(ref s_head, cell);
                Interlocked.Increment(ref s_count);
            }
        }

        // no fence here: the publication of the cell is a store that precedes the first entry of the thread in its program order, so the argument for the entry (ADR-004, lemma 1) covers it
        t_current = cell;
        return cell;
    }

    // a thread that ended holds no call (a call is a using block on its stack): the entries are 0 already, this only makes it a fact for the next owner
    private void Clear()
    {
        _first = 0;
        _nestedTop = 0;
        for (Segment? segment = _segments; segment is not null; segment = segment.Next) Array.Clear(segment.Ids);
    }
}
