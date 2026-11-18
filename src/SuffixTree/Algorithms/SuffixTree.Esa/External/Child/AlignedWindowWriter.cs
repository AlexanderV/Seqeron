using System.IO.MemoryMappedFiles;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace SuffixTree.Esa.External;

/// <summary>
/// What an <see cref="AlignedWindowWriter"/> did: windows mapped (<see cref="Opened"/>), windows mapped again after having been released (<see cref="Reopened"/>: a page that was flushed and is written again is written to storage
/// a second time, the write amplification of layout A), flushes (msync) and the bytes of the windows they covered, 32-bit writes, windows released to make room, distinct windows ever mapped.
/// </summary>
internal readonly record struct WindowStats(long Opened, long Reopened, long Flushes, long FlushedBytes, long Writes, long Evictions, long DistinctWindows);

/// <summary>
/// Writes 32-bit values at arbitrary entries of a range of a <see cref="MappedOutputFile"/> through a small set of aligned mapped windows (the child-table fields of layout A are written behind the cursor "through bounded
/// windows", ADR-004 section 3.2): a write maps the window <c>offset / W</c> (aligned, at most <c>slots</c> windows at a time, the least recently used released first); with <c>flushOnRelease</c> a released window is
/// flushed (msync), so the dirty pages of the file stay within <c>slots x W</c> plus what the kernel has not written yet. The file is created by the caller (<see cref="MappedOutputFile.Create"/>) and is not owned
/// here. Not thread-safe.
/// </summary>
internal sealed unsafe class AlignedWindowWriter : IDisposable
{
    private readonly MappedOutputFile _file;
    private readonly long _start;
    private readonly Slot[] _slots;
    private readonly int _window;
    private readonly int _shift;
    private readonly bool _flush;
    private readonly IoCounters? _counters;
    private readonly ulong[] _seen;
    private long _hotWindow = -1;
    private byte* _hotBase;
    private long _clock;
    private long _opened;
    private long _reopened;
    private long _flushes;
    private long _flushedBytes;
    private long _writes;
    private long _evictions;
    private long _distinct;
    private bool _disposed;

    private struct Slot
    {
        internal MemoryMappedViewAccessor? View;
        internal byte* Base;
        internal long Window;
        internal long Size;
        internal long LastUse;
    }

    /// <param name="file">The output file; this writer counts itself as one of its writers until it is disposed.</param>
    /// <param name="start">Offset of the range in the file.</param>
    /// <param name="length">Length of the range in bytes.</param>
    /// <param name="windowBytes">Size of a window: a power of two of at least 4.</param>
    /// <param name="slots">Windows mapped at a time.</param>
    /// <param name="flushOnRelease">msync a window when it is released.</param>
    /// <param name="counters">Receives the windows and the bytes of the writes.</param>
    internal AlignedWindowWriter(MappedOutputFile file, long start, long length, int windowBytes, int slots, bool flushOnRelease, IoCounters? counters = null)
    {
        ArgumentNullException.ThrowIfNull(file);
        ArgumentOutOfRangeException.ThrowIfNegative(start);
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        if (windowBytes < 4 || !BitOperations.IsPow2(windowBytes)) throw new ArgumentOutOfRangeException(nameof(windowBytes), windowBytes, "a window is a power of two of at least 4 bytes");
        ArgumentOutOfRangeException.ThrowIfLessThan(slots, 1);
        if (start > file.Length - length) throw new ArgumentOutOfRangeException(nameof(start), $"range [{start}, {start + length}) outside [0, {file.Length})");
        _file = file;
        _start = start;
        Length = length;
        _window = windowBytes;
        _shift = BitOperations.TrailingZeroCount(windowBytes);
        _flush = flushOnRelease;
        _counters = counters;
        _slots = new Slot[slots];
        for (int i = 0; i < slots; i++) _slots[i].Window = -1;
        long windows = (length + windowBytes - 1) >> _shift;
        _seen = new ulong[(windows + 63) / 64];
        file.WriterOpened();
    }

    /// <summary>Length of the range in bytes.</summary>
    internal long Length { get; }

    internal WindowStats Stats => new(_opened, _reopened, _flushes, _flushedBytes, _writes, _evictions, _distinct);

    /// <summary>Writes the 32-bit little-endian value at entry <paramref name="index"/> (byte offset <c>4 x index</c> of the range).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal void Write32(long index, uint value)
    {
        long offset = index << 2;
        if ((ulong)index >= (ulong)(Length >> 2)) ThrowOutside(offset, 4);
        long w = offset >> _shift;
        if (w != _hotWindow) Select(w);
        *(uint*)(_hotBase + (offset - (w << _shift))) = value;
        _writes++;
    }

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        for (int i = 0; i < _slots.Length; i++)
        {
            if (_slots[i].View is null) continue;
            if (_flush) FlushSlot(ref _slots[i]);
            Release(ref _slots[i]);
        }
        _hotWindow = -1;
        _counters?.AddMappedWrite(4 * _writes);
        _file.WriterClosed();
    }

    private void Select(long w)
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        int victim = -1;
        long oldest = long.MaxValue;
        for (int i = 0; i < _slots.Length; i++)
        {
            if (_slots[i].View is not null && _slots[i].Window == w)
            {
                _slots[i].LastUse = ++_clock;
                _hotWindow = w;
                _hotBase = _slots[i].Base;
                return;
            }
            if (_slots[i].View is null)
            {
                victim = i;
                oldest = long.MinValue;
            }
            else if (oldest != long.MinValue && _slots[i].LastUse < oldest)
            {
                victim = i;
                oldest = _slots[i].LastUse;
            }
        }
        if (_slots[victim].View is not null)
        {
            if (_flush) FlushSlot(ref _slots[victim]);
            Release(ref _slots[victim]);
            _evictions++;
        }
        Open(ref _slots[victim], w);
        _hotWindow = w;
        _hotBase = _slots[victim].Base;
    }

    private void Open(ref Slot slot, long w)
    {
        long offset = w << _shift;
        long size = Math.Min(_window, Length - offset);
        slot.View = _file.OpenView(_start + offset, size);
        byte* p = null;
        slot.View.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
        slot.Base = p + slot.View.PointerOffset;
        slot.Window = w;
        slot.Size = size;
        slot.LastUse = ++_clock;
        _opened++;
        _counters?.AddWindow(size);
        ulong bit = 1UL << (int)(w & 63);
        if ((_seen[w >> 6] & bit) != 0)
        {
            _reopened++;
        }
        else
        {
            _seen[w >> 6] |= bit;
            _distinct++;
        }
    }

    private void FlushSlot(ref Slot slot)
    {
        slot.View!.Flush();
        _flushes++;
        _flushedBytes += slot.Size;
    }

    private void Release(ref Slot slot)
    {
        slot.View!.SafeMemoryMappedViewHandle.ReleasePointer();
        slot.View.Dispose();
        slot.View = null;
        slot.Base = null;
        slot.Window = -1;
        _hotWindow = -1; // the hot window may be the one released: the next write selects again
    }

    private void ThrowOutside(long offset, long length) => throw new ArgumentOutOfRangeException(nameof(offset), $"range [{offset}, {offset + length}) outside [0, {Length})");
}
