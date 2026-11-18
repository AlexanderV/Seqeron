using System.Buffers;
using System.Runtime.InteropServices;

namespace SuffixTree.Esa.External;

/// <summary>A set of text positions, read in the order in which they were written (increasing text order for every bucket of the partition scans).</summary>
internal interface IPosSource : IDisposable
{
    /// <summary>Number of positions.</summary>
    long Count { get; }

    /// <summary>Reads positions <c>[start, start + dest.Length)</c> (the range must lie inside the set); sequential access is the fast path.</summary>
    void Read(long start, Span<uint> dest);

    /// <summary>The position at <paramref name="index"/> (sampling).</summary>
    uint At(long index);
}

/// <summary>The positions <c>0..n-1</c> without storing them (the root of the partition).</summary>
internal sealed class RangePositions(long count) : IPosSource
{
    public long Count { get; } = count;

    public void Read(long start, Span<uint> dest)
    {
        for (int i = 0; i < dest.Length; i++) dest[i] = (uint)(start + i);
    }

    public uint At(long index) => (uint)index;

    public void Dispose()
    {
        // nothing is held
    }
}

/// <summary>Positions held in an array (a bucket that never exceeded its write buffer); the array goes back to the pool when the source is disposed.</summary>
internal sealed class MemoryPositions(uint[] data, int count) : IPosSource
{
    private uint[]? _data = data;

    public long Count { get; } = count;

    public void Read(long start, Span<uint> dest) => _data.AsSpan((int)start, dest.Length).CopyTo(dest);

    public uint At(long index) => _data![index];

    public void Dispose()
    {
        uint[]? d = Interlocked.Exchange(ref _data, null);
        if (d is not null) ArrayPool<uint>.Shared.Return(d);
    }
}

/// <summary>Where the position files of one sort go: the directory, the prefix of their names, a numbering that belongs to the sort (never shared with another), the ledger and the counters.</summary>
internal sealed class PositionFiles(string directory, string prefix, ScratchLedger ledger, IoCounters? counters)
{
    private int _next;

    internal ScratchLedger Ledger { get; } = ledger;

    internal string NextPath() => Path.Combine(directory, $"{prefix}-{Interlocked.Increment(ref _next):D8}.bin");

    internal CountedFile Open(string path, FileMode mode, FileAccess access) => new(path, mode, access, counters);
}

/// <summary>Positions in a scratch file (4 bytes each, little endian); the file is deleted when the source is disposed (a scratch file lives until its consumer is done: ADR-004 Q1).</summary>
internal sealed class FilePositions(string path, long count, PositionFiles files) : IPosSource
{
    private CountedFile? _file;
    private bool _disposed;

    public long Count => count;

    public void Read(long start, Span<uint> dest)
    {
        _file ??= files.Open(path, FileMode.Open, FileAccess.Read);
        _file.ReadExactAt(MemoryMarshal.AsBytes(dest), start * 4);
    }

    public uint At(long index)
    {
        _file ??= files.Open(path, FileMode.Open, FileAccess.Read);
        Span<uint> one = stackalloc uint[1];
        _file.ReadExactAt(MemoryMarshal.AsBytes(one), index * 4);
        return one[0];
    }

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        _file?.Dispose();
        _file = null;
        files.Ledger.Remove(count * 4);
        RunStore.DeleteQuietly(path);
    }
}

/// <summary>
/// Appends positions to one bucket: it keeps up to <c>spillBytes</c> in a pooled buffer and, once the bucket outgrows the buffer, appends to a scratch file through <see cref="CountedFile"/> (write buffers: F buckets x spillBytes,
/// the "F x B_buf" term of the RAM budget in ADR-004 section 4.2). A bucket that fits stays in memory.
/// </summary>
internal sealed class PosWriter : IDisposable
{
    private readonly PositionFiles _files;
    private readonly int _capacity;
    private uint[] _buffer;
    private int _filled;
    private long _total;
    private long _written;
    private CountedFile? _file;
    private string? _path;
    private bool _completed;

    internal PosWriter(PositionFiles files, int spillBytes, int expected)
    {
        _files = files;
        _capacity = Math.Max(spillBytes / 4, 16);
        _buffer = ArrayPool<uint>.Shared.Rent(Math.Min(_capacity, Math.Max(expected, 16)));
    }

    internal long Count => _total;

    internal void Append(uint p)
    {
        if (_filled == _buffer.Length)
        {
            if (_buffer.Length < _capacity) Grow();
            else Flush();
        }
        _buffer[_filled++] = p;
        _total++;
    }

    /// <summary>Ends the bucket and hands it over as a source (the writer must not be used afterwards).</summary>
    internal IPosSource Complete()
    {
        _completed = true;
        if (_file is null)
        {
            var source = new MemoryPositions(_buffer, _filled);
            _buffer = [];
            return source;
        }
        Flush();
        _file.Dispose();
        _file = null;
        ArrayPool<uint>.Shared.Return(_buffer);
        _buffer = [];
        return new FilePositions(_path!, _total, _files);
    }

    public void Dispose()
    {
        if (_completed) return;
        _completed = true;
        _file?.Dispose();
        _file = null;
        if (_path is not null)
        {
            _files.Ledger.Remove(_written * 4);
            RunStore.DeleteQuietly(_path);
        }
        if (_buffer.Length > 0) ArrayPool<uint>.Shared.Return(_buffer);
        _buffer = [];
    }

    private void Grow()
    {
        uint[] bigger = ArrayPool<uint>.Shared.Rent(Math.Min(_capacity, _buffer.Length * 2));
        Array.Copy(_buffer, bigger, _filled);
        ArrayPool<uint>.Shared.Return(_buffer);
        _buffer = bigger;
    }

    private void Flush()
    {
        if (_file is null)
        {
            _path = _files.NextPath();
            _file = _files.Open(_path, FileMode.Create, FileAccess.ReadWrite);
        }
        if (_filled > 0)
        {
            _files.Ledger.Add((long)_filled * 4); // taken with the reservation, so that a failed write gives back exactly what was taken
            _written += _filled;
            _file.WriteAt(MemoryMarshal.AsBytes(_buffer.AsSpan(0, _filled)), (_total - _filled) * 4);
        }
        _filled = 0;
    }
}
