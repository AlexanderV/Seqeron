using Microsoft.Win32.SafeHandles;

namespace SuffixTree.Esa.External;

/// <summary>
/// A file opened for positional reads and writes (<see cref="RandomAccess"/>, i.e. pread / pwrite) that counts every byte in an optional <see cref="IoCounters"/>. The external-memory building blocks do their explicit
/// file I/O through this class, so the counters are the I/O volume of the build; accesses through memory mappings are not explicit calls and are reported separately.
/// </summary>
internal sealed class CountedFile : IDisposable
{
    private readonly SafeFileHandle _handle;
    private readonly IoCounters? _counters;

    internal CountedFile(string path, FileMode mode, FileAccess access, IoCounters? counters = null)
    {
        Path = path;
        _counters = counters;
#if ESA_BOUNDS
        if (mode is FileMode.Create or FileMode.Truncate) MappedFileRegistry.RequireNotMapped(path);
#endif
        _handle = File.OpenHandle(path, mode, access, FileShare.ReadWrite, FileOptions.None);
    }

    internal string Path { get; }

    internal long Length => RandomAccess.GetLength(_handle);

    /// <summary>Reads up to <c>buffer.Length</c> bytes at <paramref name="offset"/> (fewer only at the end of the file); returns the count.</summary>
    internal int ReadAt(Span<byte> buffer, long offset)
    {
        int total = 0;
        while (total < buffer.Length)
        {
            int got = RandomAccess.Read(_handle, buffer[total..], offset + total);
            if (got <= 0) break;
            total += got;
            _counters?.AddRead(got);
        }
        return total;
    }

    /// <summary>Reads exactly <c>buffer.Length</c> bytes or throws <see cref="EndOfStreamException"/>.</summary>
    internal void ReadExactAt(Span<byte> buffer, long offset)
    {
        if (ReadAt(buffer, offset) != buffer.Length) throw new EndOfStreamException($"{Path}: short read at offset {offset}");
    }

    /// <summary>Writes all of <paramref name="buffer"/> at <paramref name="offset"/>.</summary>
    internal void WriteAt(ReadOnlySpan<byte> buffer, long offset)
    {
        RandomAccess.Write(_handle, buffer, offset);
        _counters?.AddWrite(buffer.Length);
    }

    internal void SetLength(long length) => RandomAccess.SetLength(_handle, length);

    /// <summary>Forces the written data to storage (fsync): bounds the dirty pages of a scratch file that is written once and read back later.</summary>
    internal void Flush()
    {
        RandomAccess.FlushToDisk(_handle);
        _counters?.AddFlush();
    }

    public void Dispose() => _handle.Dispose();
}
