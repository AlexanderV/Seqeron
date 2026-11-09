using System.Buffers;
using System.Runtime.InteropServices;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Build;

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
    }
}

/// <summary>Positions held in a managed array (memory mode, or a bucket that never exceeded its write buffer).</summary>
internal sealed class MemoryPositions(uint[] data, int count, bool pooled) : IPosSource
{
    private uint[]? _data = data;

    public long Count { get; } = count;

    public void Read(long start, Span<uint> dest) => _data.AsSpan((int)start, dest.Length).CopyTo(dest);

    public uint At(long index) => _data![index];

    public void Dispose()
    {
        uint[]? d = Interlocked.Exchange(ref _data, null);
        if (d is not null && pooled) ArrayPool<uint>.Shared.Return(d);
    }
}

/// <summary>Positions in a scratch file (4 bytes each, little endian); the file is deleted when the source is disposed (a scratch file lives until its consumer is done: ADR-004 Q1).</summary>
internal sealed class FilePositions(string path, long count) : IPosSource
{
    private readonly long _count = count;
    private CountedFile? _file;
    private bool _disposed;

    public long Count => _count;

    public void Read(long start, Span<uint> dest)
    {
        _file ??= new CountedFile(path, FileMode.Open, FileAccess.Read);
        _file.ReadExactAt(MemoryMarshal.AsBytes(dest), start * 4);
    }

    public uint At(long index)
    {
        _file ??= new CountedFile(path, FileMode.Open, FileAccess.Read);
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
        ScratchTracker.Remove(_count * 4);
        try
        {
            File.Delete(path);
        }
        catch (IOException)
        {
            // best effort: the scratch directory is removed by the owner at the end of the build
        }
    }
}

/// <summary>
/// Appends positions to one bucket. Memory mode (<c>dir</c> null) keeps everything in an array. File mode buffers <c>spillBytes</c> and, once the bucket outgrows the buffer, appends to a scratch file
/// through <see cref="CountedFile"/> (write buffers: F buckets x spillBytes, the "F x B_buf" term of the RAM budget in ADR-004 section 4.2).
/// </summary>
internal sealed class PosWriter : IDisposable
{
    private readonly string? _dir;
    private readonly int _capacity;
    private uint[] _buffer;
    private int _filled;
    private long _total;
    private CountedFile? _file;
    private string? _path;
    private bool _completed;

    internal PosWriter(string? dir, int spillBytes, int expected)
    {
        _dir = dir;
        _capacity = dir is null ? int.MaxValue : Math.Max(spillBytes / 4, 16);
        _buffer = ArrayPool<uint>.Shared.Rent(dir is null ? Math.Max(expected, 16) : Math.Min(_capacity, Math.Max(expected, 16)));
    }

    internal long Count => _total;

    internal void Append(uint p)
    {
        if (_filled == _buffer.Length)
        {
            if (_dir is null || _buffer.Length < _capacity) Grow();
            else Flush();
        }
        _buffer[_filled++] = p;
        _total++;
    }

    internal IPosSource Complete()
    {
        _completed = true;
        if (_file is null)
        {
            var source = new MemoryPositions(_buffer, _filled, pooled: true);
            _buffer = [];
            return source;
        }
        Flush();
        _file.Dispose();
        _file = null;
        ArrayPool<uint>.Shared.Return(_buffer);
        _buffer = [];
        return new FilePositions(_path!, _total);
    }

    public void Dispose()
    {
        if (_completed) return;
        _file?.Dispose();
        if (_path is not null)
        {
            ScratchTracker.Remove((_total - _filled) * 4);
            try
            {
                File.Delete(_path);
            }
            catch (IOException)
            {
                // ignored: see FilePositions.Dispose
            }
        }
        if (_buffer.Length > 0) ArrayPool<uint>.Shared.Return(_buffer);
        _buffer = [];
    }

    private void Grow()
    {
        int target = _dir is null ? _buffer.Length * 2 : Math.Min(_capacity, _buffer.Length * 2);
        uint[] bigger = ArrayPool<uint>.Shared.Rent(target);
        Array.Copy(_buffer, bigger, _filled);
        ArrayPool<uint>.Shared.Return(_buffer);
        _buffer = bigger;
    }

    private void Flush()
    {
        if (_file is null)
        {
            _path = Path.Combine(_dir!, $"pos-{Interlocked.Increment(ref s_fileCounter):D8}.bin");
            _file = new CountedFile(_path, FileMode.Create, FileAccess.ReadWrite);
        }
        if (_filled > 0)
        {
            _file.WriteAt(MemoryMarshal.AsBytes(_buffer.AsSpan(0, _filled)), (_total - _filled) * 4);
            ScratchTracker.Add((long)_filled * 4);
        }
        _filled = 0;
    }

    private static long s_fileCounter;
}
