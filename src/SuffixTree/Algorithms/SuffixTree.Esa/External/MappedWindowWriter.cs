#pragma warning disable S2325, CA1822, S1144, IDE0051, IDE0060, S1172 // red stub of work package 1b.7a: every member throws NotImplementedException; the implementation replaces this file

namespace SuffixTree.Esa.External;

/// <summary>Writes a file of known length through bounded sliding windows of a memory mapping. Not thread-safe.</summary>
internal sealed unsafe class MappedWindowWriter : IDisposable
{
    /// <summary>Default window: 1 MiB (the .NET view accessor flushes every window when it is disposed, so small windows multiply the device writes).</summary>
    internal const int DefaultWindowBytes = 1 << 20;

    internal MappedWindowWriter(string path, long length, int windowBytes = DefaultWindowBytes, IoCounters? counters = null) => throw new NotImplementedException();

    internal string Path => throw new NotImplementedException();

    internal long Length => throw new NotImplementedException();

    /// <summary>Offset of the next <see cref="Append"/>.</summary>
    internal long Position => throw new NotImplementedException();

    internal void WriteAt(long offset, ReadOnlySpan<byte> data) => throw new NotImplementedException();

    internal void Append(ReadOnlySpan<byte> data) => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
