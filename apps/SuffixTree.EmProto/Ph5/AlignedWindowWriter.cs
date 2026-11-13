using System.IO.MemoryMappedFiles;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace SuffixTree.EmProto.Ph5;

/// <summary>
/// What an <see cref="AlignedWindowWriter"/> did: windows mapped (<see cref="Opened"/>), windows mapped again after having been released (<see cref="Reopened"/>: a page that was flushed and is written again is written to storage
/// a second time, the write amplification of layout A), flushes (msync) and the bytes of the windows they covered, 32-bit writes, windows released to make room, distinct windows ever mapped.
/// </summary>
internal readonly record struct WindowStats(long Opened, long Reopened, long Flushes, long FlushedBytes, long Writes, long Evictions, long DistinctWindows);

/// <summary>
/// Writes a file of known length through a small set of aligned mapped windows (the categorical rule of the owner, 2026-11-04: stream into memory-mapped files, no "build in RAM and then write"; ADR-004 section 3.2: the
/// child-table fields of layout A are written behind the cursor "through bounded windows"). The file is created at its final length; a write maps the window <c>offset / W</c> (aligned, at most <c>slots</c> windows at a
/// time, the least recently used released first); with <c>flushOnRelease</c> a released window is flushed (<see cref="System.IO.MemoryMappedFiles.MemoryMappedViewAccessor.Flush"/>, msync), so the dirty pages of the
/// file stay within <c>slots x W</c> plus what the kernel has not written yet. The same class writes layout A (random writes behind the cursor), the PSV (forward) and the NSV (backward). Not thread-safe.
/// </summary>
internal sealed unsafe class AlignedWindowWriter : IDisposable
{
    private readonly MemoryMappedFile? _file;
    private readonly Slot[] _slots;
    private readonly int _window;
    private readonly int _shift;
    private readonly bool _flush;
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

    internal AlignedWindowWriter(string path, long length, int windowBytes, int slots, bool flushOnRelease)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        if (windowBytes < 4 || !BitOperations.IsPow2(windowBytes)) throw new ArgumentOutOfRangeException(nameof(windowBytes), windowBytes, "a window is a power of two of at least 4 bytes");
        ArgumentOutOfRangeException.ThrowIfLessThan(slots, 1);
        Length = length;
        _window = windowBytes;
        _shift = BitOperations.TrailingZeroCount(windowBytes);
        _flush = flushOnRelease;
        _slots = new Slot[slots];
        for (int i = 0; i < slots; i++) _slots[i].Window = -1;
        long windows = (length + windowBytes - 1) >> _shift;
        _seen = new ulong[(windows + 63) / 64];
        using (var fs = new FileStream(path, FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite)) fs.SetLength(length);
        if (length > 0) _file = MemoryMappedFile.CreateFromFile(path, FileMode.Open, null, 0, MemoryMappedFileAccess.ReadWrite);
    }

    /// <summary>Length of the file in bytes.</summary>
    internal long Length { get; }

    internal WindowStats Stats => new(_opened, _reopened, _flushes, _flushedBytes, _writes, _evictions, _distinct);

    /// <summary>Writes the 32-bit little-endian value at entry <paramref name="index"/> (byte offset <c>4 x index</c>).</summary>
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

    /// <summary>Writes <paramref name="data"/> at <paramref name="offset"/> (the range may span windows).</summary>
    internal void Write(long offset, ReadOnlySpan<byte> data)
    {
        if (offset < 0 || offset + data.Length > Length) ThrowOutside(offset, data.Length);
        while (!data.IsEmpty)
        {
            long w = offset >> _shift;
            if (w != _hotWindow) Select(w);
            long inside = offset - (w << _shift);
            int take = (int)Math.Min(data.Length, _slots[HotSlot()].Size - inside);
            data[..take].CopyTo(new Span<byte>(_hotBase + inside, take));
            offset += take;
            data = data[take..];
        }
    }

    /// <summary>Flushes every mapped window (they stay mapped).</summary>
    internal void Flush()
    {
        for (int i = 0; i < _slots.Length; i++)
        {
            if (_slots[i].View is not null) FlushSlot(ref _slots[i]);
        }
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
        _file?.Dispose();
    }

    private int HotSlot()
    {
        for (int i = 0; i < _slots.Length; i++)
        {
            if (_slots[i].Window == _hotWindow && _slots[i].View is not null) return i;
        }
        throw new InvalidOperationException("no hot window");
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
        long start = w << _shift;
        long size = Math.Min(_window, Length - start);
        slot.View = _file!.CreateViewAccessor(start, size, MemoryMappedFileAccess.Write);
        byte* p = null;
        slot.View.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
        slot.Base = p + slot.View.PointerOffset;
        slot.Window = w;
        slot.Size = size;
        slot.LastUse = ++_clock;
        _opened++;
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
