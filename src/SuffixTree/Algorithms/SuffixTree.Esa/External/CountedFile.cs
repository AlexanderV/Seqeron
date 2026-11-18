#pragma warning disable S2325, CA1822, S1144, IDE0051, IDE0060, S1172 // red stub of work package 1b.7a: every member throws NotImplementedException; the implementation replaces this file

namespace SuffixTree.Esa.External;

/// <summary>A file opened for positional reads and writes (pread / pwrite) that counts every byte in an optional <see cref="IoCounters"/>.</summary>
internal sealed class CountedFile : IDisposable
{
    internal CountedFile(string path, FileMode mode, FileAccess access, IoCounters? counters = null) => throw new NotImplementedException();

    internal string Path => throw new NotImplementedException();

    internal long Length => throw new NotImplementedException();

    /// <summary>Reads up to <c>buffer.Length</c> bytes at <paramref name="offset"/> (fewer only at the end of the file); returns the count.</summary>
    internal int ReadAt(Span<byte> buffer, long offset) => throw new NotImplementedException();

    /// <summary>Reads exactly <c>buffer.Length</c> bytes or throws <see cref="EndOfStreamException"/>.</summary>
    internal void ReadExactAt(Span<byte> buffer, long offset) => throw new NotImplementedException();

    internal void WriteAt(ReadOnlySpan<byte> buffer, long offset) => throw new NotImplementedException();

    internal void SetLength(long length) => throw new NotImplementedException();

    /// <summary>Forces the written data to storage (fsync).</summary>
    internal void Flush() => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
