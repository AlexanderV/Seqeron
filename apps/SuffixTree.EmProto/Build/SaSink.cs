using System.IO.MemoryMappedFiles;
using System.Security.Cryptography;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Build;

/// <summary>Destination of the suffix array (n + 1 entries of 4 bytes, <c>SA[0] = n</c>): every bucket writes its slice at its final offset, exactly once, from any thread.</summary>
internal interface ISaSink : IDisposable
{
    /// <summary>Writes <paramref name="data"/> at entry index <paramref name="offset"/>; slices of different calls never overlap.</summary>
    void Write(long offset, ReadOnlySpan<uint> data);
}

/// <summary>The suffix array in native memory (E1, E2 and the tests).</summary>
internal sealed unsafe class MemorySaSink : ISaSink
{
    private readonly NativeBlock _block;

    internal MemorySaSink(long entries)
    {
        Entries = entries;
        _block = new NativeBlock(entries * 4);
    }

    internal long Entries { get; }

    internal uint* Data => (uint*)_block.Pointer;

    public void Write(long offset, ReadOnlySpan<uint> data)
    {
        if (offset < 0 || offset + data.Length > Entries) throw new ArgumentOutOfRangeException(nameof(offset), $"slice [{offset}, {offset + data.Length}) outside [0, {Entries})");
        data.CopyTo(new Span<uint>(Data + offset, data.Length));
    }

    internal uint[] ToArray() => new ReadOnlySpan<uint>(Data, checked((int)Entries)).ToArray();

    /// <summary>SHA-256 (lower-case hex) of entries <c>1..n</c> as little-endian int32: the libsais hash of <c>analysis/sa_vectors.tsv</c> (this machine is little endian).</summary>
    internal string HashLikeLibsais() => SaHash.Compute(Data + 1, Entries - 1);

    public void Dispose() => _block.Dispose();
}

/// <summary>The suffix array in a file, written through bounded memory-mapped windows: one view per slice (at most the leaf capacity x 4 bytes), disposed right after the copy, never flushed explicitly.</summary>
internal sealed unsafe class MappedSaSink : ISaSink
{
    private readonly MemoryMappedFile _file;
    private readonly long _entries;

    internal MappedSaSink(string path, long entries)
    {
        _entries = entries;
        Path = path;
        using (var fs = new FileStream(path, FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite)) fs.SetLength(entries * 4);
        _file = MemoryMappedFile.CreateFromFile(path, FileMode.Open, null, 0, MemoryMappedFileAccess.ReadWrite);
    }

    internal string Path { get; }

    /// <summary>Largest window created so far, in bytes (reported as the bound of the mapped write window).</summary>
    internal long MaxWindowBytes => Interlocked.Read(ref _maxWindow);

    private long _maxWindow;

    public void Write(long offset, ReadOnlySpan<uint> data)
    {
        if (data.IsEmpty) return;
        if (offset < 0 || offset + data.Length > _entries) throw new ArgumentOutOfRangeException(nameof(offset), $"slice [{offset}, {offset + data.Length}) outside [0, {_entries})");
        long bytes = (long)data.Length * 4;
        Harness.Atomic.RaiseTo(ref _maxWindow, bytes);
        using MemoryMappedViewAccessor view = _file.CreateViewAccessor(offset * 4, bytes, MemoryMappedFileAccess.Write);
        byte* p = null;
        view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
        try
        {
            data.CopyTo(new Span<uint>(p + view.PointerOffset, data.Length));
        }
        finally
        {
            view.SafeMemoryMappedViewHandle.ReleasePointer();
        }
    }

    public void Dispose() => _file.Dispose();
}

/// <summary>The libsais-compatible hash of a suffix array.</summary>
internal static unsafe class SaHash
{
    /// <summary>SHA-256 (lower-case hex) of <paramref name="count"/> little-endian int32 values at <paramref name="p"/>.</summary>
    internal static string Compute(uint* p, long count)
    {
        using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        const int chunk = 1 << 20;
        for (long done = 0; done < count; done += chunk)
        {
            int take = (int)Math.Min(chunk, count - done);
            hash.AppendData(new ReadOnlySpan<byte>(p + done, take * 4));
        }
        return Convert.ToHexStringLower(hash.GetHashAndReset());
    }

    /// <summary>The same hash over a file written by <see cref="MappedSaSink"/>: entries <c>1..n</c> (the first four bytes, SA[0], are skipped).</summary>
    internal static string ComputeFile(string path)
    {
        using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20, FileOptions.SequentialScan);
        fs.Seek(4, SeekOrigin.Begin);
        var buf = new byte[1 << 20];
        int got;
        while ((got = fs.Read(buf, 0, buf.Length)) > 0) hash.AppendData(buf.AsSpan(0, got));
        return Convert.ToHexStringLower(hash.GetHashAndReset());
    }
}
