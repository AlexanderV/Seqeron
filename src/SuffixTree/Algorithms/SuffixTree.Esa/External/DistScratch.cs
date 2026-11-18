#pragma warning disable S2325, CA1822, S1144, IDE0051, IDE0060, S1172 // red stub of work package 1b.7a: every member throws NotImplementedException; the implementation replaces this file

namespace SuffixTree.Esa.External;

/// <summary>The scratch side of the DIST primitive: P chunk files of fixed-size tuples.</summary>
internal sealed unsafe class DistScratch<T> : IDisposable
    where T : unmanaged
{
    internal DistScratch(string dir, string tag, long chunks, IoCounters? counters = null) => throw new NotImplementedException();

    internal long Chunks => throw new NotImplementedException();

    internal long Entries(long chunk) => throw new NotImplementedException();

    internal string PathOf(long chunk) => throw new NotImplementedException();

    internal void Append(long chunk, ReadOnlySpan<T> items) => throw new NotImplementedException();

    internal DistReader<T> OpenReader(long chunk) => throw new NotImplementedException();

    internal void Delete(long chunk) => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}

/// <summary>Sequential reader of one chunk of a <see cref="DistScratch{T}"/>.</summary>
internal sealed unsafe class DistReader<T> : IDisposable
    where T : unmanaged
{
    internal long Count => throw new NotImplementedException();

    internal int Read(Span<T> buffer) => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}

/// <summary>One task's write side of a distribution: a buffer per chunk, flushed when full.</summary>
internal sealed unsafe class DistWriter<T>
    where T : unmanaged
{
    internal DistWriter(DistScratch<T> scratch, BudgetedArena arena, long firstChunk, long chunkCount, int bufferBytes) => throw new NotImplementedException();

    internal int Capacity => throw new NotImplementedException();

    internal void Add(long chunk, in T item) => throw new NotImplementedException();

    /// <summary>Entries of <paramref name="chunk"/> that are still in the buffer.</summary>
    internal int Pending(long chunk) => throw new NotImplementedException();

    /// <summary>The buffered entries of <paramref name="chunk"/> (valid until the next <see cref="Add"/> to it).</summary>
    internal ReadOnlySpan<T> PendingSpan(long chunk) => throw new NotImplementedException();

    internal void FlushAll() => throw new NotImplementedException();
}
