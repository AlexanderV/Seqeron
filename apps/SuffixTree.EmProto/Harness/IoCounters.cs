using Microsoft.Win32.SafeHandles;

namespace SuffixTree.EmProto.Harness;

/// <summary>Process-wide byte and call counters of every explicit file read and write of the prototype (E0: they must match <c>/proc/self/io</c> within 1 %).</summary>
internal static class IoCounters
{
    private static long _readBytes;
    private static long _writeBytes;
    private static long _readCalls;
    private static long _writeCalls;

    internal static long ReadBytes => Interlocked.Read(ref _readBytes);

    internal static long WriteBytes => Interlocked.Read(ref _writeBytes);

    internal static long ReadCalls => Interlocked.Read(ref _readCalls);

    internal static long WriteCalls => Interlocked.Read(ref _writeCalls);

    internal static void Reset()
    {
        Interlocked.Exchange(ref _readBytes, 0);
        Interlocked.Exchange(ref _writeBytes, 0);
        Interlocked.Exchange(ref _readCalls, 0);
        Interlocked.Exchange(ref _writeCalls, 0);
    }

    internal static void AddRead(long bytes)
    {
        Interlocked.Add(ref _readBytes, bytes);
        Interlocked.Increment(ref _readCalls);
    }

    internal static void AddWrite(long bytes)
    {
        Interlocked.Add(ref _writeBytes, bytes);
        Interlocked.Increment(ref _writeCalls);
    }
}

/// <summary>
/// A file opened for positional reads and writes (<see cref="RandomAccess"/>, i.e. pread / pwrite) that counts every byte in <see cref="IoCounters"/>.
/// The whole prototype does its explicit file I/O through this class, so the counters are the I/O volume of the build; accesses through memory
/// mappings are not explicit calls and are reported separately (page-cache and block-layer numbers come from <c>/proc/self/io</c> and <c>/proc/meminfo</c>).
/// </summary>
internal sealed class CountedFile : IDisposable
{
    private readonly SafeFileHandle _handle;

    internal CountedFile(string path, FileMode mode, FileAccess access)
    {
        Path = path;
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
            IoCounters.AddRead(got);
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
        IoCounters.AddWrite(buffer.Length);
    }

    internal void SetLength(long length) => RandomAccess.SetLength(_handle, length);

    /// <summary>Forces the written data to storage (fsync): bounds the dirty pages of a scratch file that is written once and read back later.</summary>
    internal void Flush() => RandomAccess.FlushToDisk(_handle);

    public void Dispose() => _handle.Dispose();
}
