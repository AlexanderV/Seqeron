using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Ph4;

/// <summary>Distribution tuple of the base mode, written by B1: the position <c>SA[i]</c>, its predecessor <c>SA[i-1]</c> (the Phi value) and the rank <c>i</c>; 12 bytes.</summary>
[StructLayout(LayoutKind.Sequential)]
internal readonly struct PhiTuple(uint position, uint phi, uint rank)
{
    internal readonly uint Position = position;
    internal readonly uint Phi = phi;
    internal readonly uint Rank = rank;
}

/// <summary>Distribution tuple of the bounded mode (the rank is not needed, the gather uses the suffix array again): position and Phi value; 8 bytes.</summary>
[StructLayout(LayoutKind.Sequential)]
internal readonly struct PhiPair(uint position, uint phi)
{
    internal readonly uint Position = position;
    internal readonly uint Phi = phi;
}

/// <summary>Distribution tuple of the second pass of the base mode, written by B2: a rank and the LCP value at that rank (<c>PLCP[SA[i]]</c>); 8 bytes.</summary>
[StructLayout(LayoutKind.Sequential)]
internal readonly struct RankTuple(uint rank, uint value)
{
    internal readonly uint Rank = rank;
    internal readonly uint Value = value;
}

/// <summary>
/// The scratch side of the DIST primitive of ADR-004 section 3.3, generalised from <see cref="Dcs.SaInverter"/> for the tuple type: P chunk files <c>tag-chunk.bin</c> of fixed-size tuples, each appended
/// to by one or several <see cref="SpoolWriter{T}"/> (the append offset of a chunk is reserved atomically, so tasks never overlap), read back chunk by chunk in order through <see cref="SpoolReader{T}"/> and
/// deleted as soon as the chunk has been consumed (ADR-004 Q1: scratch is released stage by stage). Every explicit read and write goes through <see cref="CountedFile"/>, every byte on disk through
/// <see cref="ScratchTracker"/>. Disposing removes the files that are left (a failed build).
/// </summary>
internal sealed unsafe class Spool<T> : IDisposable
    where T : unmanaged
{
    private readonly string _dir;
    private readonly string _tag;
    private readonly long[] _entries;

    internal Spool(string dir, string tag, long chunks)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(chunks);
        _dir = dir;
        _tag = tag;
        _entries = new long[chunks];
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
        using (var f = new CountedFile(PathOf(chunk), FileMode.OpenOrCreate, FileAccess.Write))
        {
            f.WriteAt(MemoryMarshal.AsBytes(items), at * sizeof(T));
        }
        ScratchTracker.Add((long)items.Length * sizeof(T));
    }

    /// <summary>Opens chunk <paramref name="chunk"/> for sequential reading.</summary>
    internal SpoolReader<T> OpenReader(long chunk)
    {
        long count = Entries(chunk);
        return new SpoolReader<T>(count > 0 ? new CountedFile(PathOf(chunk), FileMode.Open, FileAccess.Read) : null, count);
    }

    /// <summary>Removes the file of <paramref name="chunk"/> and gives its bytes back to the scratch tracker.</summary>
    internal void Delete(long chunk)
    {
        long count = Interlocked.Exchange(ref _entries[chunk], 0);
        if (count == 0) return;
        ScratchTracker.Remove(count * sizeof(T));
        File.Delete(PathOf(chunk));
    }

    public void Dispose()
    {
        for (long c = 0; c < _entries.Length; c++) Delete(c);
    }
}

/// <summary>Sequential reader of one chunk of a <see cref="Spool{T}"/>.</summary>
internal sealed unsafe class SpoolReader<T> : IDisposable
    where T : unmanaged
{
    private readonly CountedFile? _file;
    private long _position;

    internal SpoolReader(CountedFile? file, long count)
    {
        _file = file;
        Count = count;
    }

    /// <summary>Entries in the chunk.</summary>
    internal long Count { get; }

    /// <summary>Reads the next entries into <paramref name="buffer"/>; returns the number read (0 at the end).</summary>
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
/// One task's write side of a distribution: a buffer of at most <c>bufferBytes</c> per chunk (all taken from one arena block, so the sum of the P buffers is charged to the budget), flushed to the chunk
/// file of the <see cref="Spool{T}"/> when it is full and by <see cref="FlushAll"/>. Not thread-safe: every task has its own writer, the files are shared.
/// </summary>
internal sealed unsafe class SpoolWriter<T>
    where T : unmanaged
{
    private readonly Spool<T> _spool;
    private readonly Buf<T> _buffers;
    private readonly int[] _filled;
    private readonly long _firstChunk;
    private readonly int _capacity;

    /// <param name="spool">The files.</param>
    /// <param name="arena">Owner of the buffer block.</param>
    /// <param name="firstChunk">The first chunk this writer will see.</param>
    /// <param name="chunkCount">Number of consecutive chunks from <paramref name="firstChunk"/>.</param>
    /// <param name="bufferBytes">Buffer per chunk (at least one tuple).</param>
    internal SpoolWriter(Spool<T> spool, BudgetedArena arena, long firstChunk, long chunkCount, int bufferBytes)
    {
        _spool = spool;
        _firstChunk = firstChunk;
        _capacity = Math.Max(1, bufferBytes / sizeof(T));
        _filled = new int[chunkCount];
        _buffers = arena.Allocate<T>(checked((int)(chunkCount * _capacity)), zero: false, "spool.buffers");
    }

    /// <summary>Tuples a chunk buffer holds before it is flushed.</summary>
    internal int Capacity => _capacity;

    /// <summary>Appends one tuple to the buffer of <paramref name="chunk"/>.</summary>
    internal void Add(long chunk, in T item)
    {
        int local = (int)(chunk - _firstChunk);
        int filled = _filled[local];
        _buffers[(local * _capacity) + filled] = item;
        if (++filled == _capacity)
        {
            Flush(local, filled);
            filled = 0;
        }
        _filled[local] = filled;
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

    private void Flush(int local, int count) => _spool.Append(_firstChunk + local, _buffers.Slice(local * _capacity, count).AsSpan(count));
}
