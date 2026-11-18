using System.Runtime.InteropServices;

namespace SuffixTree.Esa.External;

/// <summary>
/// The scratch side of the DIST primitive of ADR-004 section 3.3, generalised in the tuple type (the prototype's <c>Spool</c>): P chunk files <c>tag-chunk.bin</c> of fixed-size tuples, each appended to by one or
/// several <see cref="DistWriter{T}"/> (the append offset of a chunk is reserved atomically, so tasks never overlap), read back chunk by chunk in order through <see cref="DistReader{T}"/> and deleted as soon as the
/// chunk has been consumed (ADR-004 Q1: scratch is released stage by stage). Every explicit read and write goes through <see cref="CountedFile"/>, every byte on disk through the <see cref="IoCounters"/> scratch
/// level. Disposing removes the files that are left (a failed build).
/// </summary>
internal sealed unsafe class DistScratch<T> : IDisposable
    where T : unmanaged
{
    private readonly string _dir;
    private readonly string _tag;
    private readonly long[] _entries;
    private readonly IoCounters? _counters;

    internal DistScratch(string dir, string tag, long chunks, IoCounters? counters = null)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(chunks);
        _dir = dir;
        _tag = tag;
        _entries = new long[chunks];
        _counters = counters;
    }

    /// <summary>Number of chunk files (some may never be created).</summary>
    internal long Chunks => _entries.Length;

    /// <summary>Entries written to chunk <paramref name="chunk"/> so far.</summary>
    internal long Entries(long chunk) => Interlocked.Read(ref _entries[chunk]);

    internal string PathOf(long chunk) => Path.Combine(_dir, $"{_tag}-{chunk:D8}.bin");

    /// <summary>Appends <paramref name="items"/> to the file of <paramref name="chunk"/> (thread-safe; the file is created by the first append).</summary>
    internal void Append(long chunk, ReadOnlySpan<T> items)
    {
        if (items.IsEmpty) return;
        long at = Interlocked.Add(ref _entries[chunk], items.Length) - items.Length;
        _counters?.AddScratch((long)items.Length * sizeof(T)); // accounted with the reservation, so that Delete gives back exactly what was taken, also after a failed write
        using var f = new CountedFile(PathOf(chunk), FileMode.OpenOrCreate, FileAccess.Write, _counters);
        f.WriteAt(MemoryMarshal.AsBytes(items), at * sizeof(T));
    }

    /// <summary>Opens chunk <paramref name="chunk"/> for sequential reading.</summary>
    internal DistReader<T> OpenReader(long chunk)
    {
        long count = Entries(chunk);
        return new DistReader<T>(count > 0 ? new CountedFile(PathOf(chunk), FileMode.Open, FileAccess.Read, _counters) : null, count);
    }

    /// <summary>Removes the file of <paramref name="chunk"/> and gives its bytes back to the scratch accounting (once).</summary>
    internal void Delete(long chunk)
    {
        long count = Interlocked.Exchange(ref _entries[chunk], 0);
        if (count == 0) return;
        _counters?.RemoveScratch(count * sizeof(T));
        File.Delete(PathOf(chunk));
    }

    public void Dispose()
    {
        for (long c = 0; c < _entries.Length; c++) Delete(c);
    }
}

/// <summary>Sequential reader of one chunk of a <see cref="DistScratch{T}"/>.</summary>
internal sealed unsafe class DistReader<T> : IDisposable
    where T : unmanaged
{
    private readonly CountedFile? _file;
    private long _position;

    internal DistReader(CountedFile? file, long count)
    {
        _file = file;
        Count = count;
    }

    /// <summary>Entries in the chunk.</summary>
    internal long Count { get; }

    /// <summary>Reads the next entries into <paramref name="buffer"/>; returns the number read (0 at the end). A chunk file that is shorter than its count throws <see cref="EndOfStreamException"/>.</summary>
    internal int Read(Span<T> buffer)
    {
        int take = (int)Math.Min(buffer.Length, Count - _position);
        if (take <= 0) return 0;
        _file!.ReadExactAt(MemoryMarshal.AsBytes(buffer[..take]), _position * sizeof(T));
        _position += take;
        return take;
    }

    public void Dispose() => _file?.Dispose();
}

/// <summary>
/// One task's write side of a distribution (the prototype's <c>SpoolWriter</c>): a buffer of at most <c>bufferBytes</c> per chunk (all taken from one arena block, so the sum of the P buffers is charged to the
/// budget), flushed to the chunk file of the <see cref="DistScratch{T}"/> when it is full and by <see cref="FlushAll"/>. What is still in a buffer is readable through <see cref="PendingSpan"/>, so a chunk that
/// never outgrows its buffer need not touch the disk. Not thread-safe: every task has its own writer, the files are shared.
/// </summary>
internal sealed unsafe class DistWriter<T>
    where T : unmanaged
{
    private readonly DistScratch<T> _scratch;
    private readonly Buf<T> _buffers;
    private readonly int[] _filled;
    private readonly long _firstChunk;
    private readonly int _capacity;

    /// <param name="scratch">The files.</param>
    /// <param name="arena">Owner of the buffer block.</param>
    /// <param name="firstChunk">The first chunk this writer will see.</param>
    /// <param name="chunkCount">Number of consecutive chunks from <paramref name="firstChunk"/>.</param>
    /// <param name="bufferBytes">Buffer per chunk (at least one tuple).</param>
    internal DistWriter(DistScratch<T> scratch, BudgetedArena arena, long firstChunk, long chunkCount, int bufferBytes)
    {
        ArgumentNullException.ThrowIfNull(scratch);
        ArgumentNullException.ThrowIfNull(arena);
        ArgumentOutOfRangeException.ThrowIfNegative(firstChunk);
        ArgumentOutOfRangeException.ThrowIfNegative(chunkCount);
        _scratch = scratch;
        _firstChunk = firstChunk;
        _capacity = Math.Max(1, bufferBytes / sizeof(T));
        _filled = new int[chunkCount];
        _buffers = arena.Allocate<T>(checked((int)(chunkCount * _capacity)), zero: false, "dist.tuples");
    }

    /// <summary>Tuples a chunk buffer holds before it is flushed.</summary>
    internal int Capacity => _capacity;

    /// <summary>Appends one tuple to the buffer of <paramref name="chunk"/>.</summary>
    internal void Add(long chunk, in T item)
    {
        int local = Local(chunk);
        int filled = _filled[local];
        _buffers[(local * _capacity) + filled] = item;
        if (++filled == _capacity)
        {
            Flush(local, filled);
            filled = 0;
        }
        _filled[local] = filled;
    }

    /// <summary>Entries of <paramref name="chunk"/> that are still in the buffer.</summary>
    internal int Pending(long chunk) => _filled[Local(chunk)];

    /// <summary>The buffered entries of <paramref name="chunk"/> (valid until the next <see cref="Add"/> to it).</summary>
    internal ReadOnlySpan<T> PendingSpan(long chunk)
    {
        int local = Local(chunk);
        int filled = _filled[local];
        return _buffers.Slice(local * _capacity, filled).AsSpan(filled);
    }

    /// <summary>Writes out every partly filled buffer.</summary>
    internal void FlushAll()
    {
        for (int local = 0; local < _filled.Length; local++)
        {
            if (_filled[local] == 0) continue;
            Flush(local, _filled[local]);
            _filled[local] = 0;
        }
    }

    private int Local(long chunk)
    {
        long local = chunk - _firstChunk;
        if ((ulong)local >= (ulong)_filled.Length) throw new ArgumentOutOfRangeException(nameof(chunk), $"chunk {chunk} is outside [{_firstChunk}, {_firstChunk + _filled.Length})");
        return (int)local;
    }

    private void Flush(int local, int count) => _scratch.Append(_firstChunk + local, _buffers.Slice(local * _capacity, count).AsSpan(count));
}
