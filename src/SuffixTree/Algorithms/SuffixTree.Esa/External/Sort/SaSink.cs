using System.Runtime.InteropServices;

namespace SuffixTree.Esa.External;

/// <summary>Where the leaf buckets write their slices of the suffix array (entry offsets, 4 bytes per entry): every slice at its final offset, exactly once, from any worker.</summary>
internal interface ISaSliceSink
{
    /// <summary>A writer for one worker (not thread-safe; the worker disposes it).</summary>
    ISaSliceWriter CreateWriter();
}

/// <summary>The writer of one worker of an <see cref="ISaSliceSink"/>.</summary>
internal interface ISaSliceWriter : IDisposable
{
    /// <summary>Writes <paramref name="data"/> at the entry <paramref name="entryOffset"/>; slices of different calls never overlap.</summary>
    void Write(long entryOffset, ReadOnlySpan<uint> data);
}

/// <summary>
/// The suffix array file: a <see cref="MappedOutputFile"/> of <c>4 (n + 1)</c> bytes, created at its final length (so a full disk is an <see cref="IOException"/> at creation), and one sliding-window writer per worker over the
/// whole range (disjoint slices are written at their offsets; the window is at least 1 MiB in production, the E9 rule, and follows the worker from slice to slice).
/// </summary>
internal sealed class MappedSaSliceSink(MappedOutputFile file, int windowBytes, IoCounters? counters) : ISaSliceSink
{
    public ISaSliceWriter CreateWriter() => new Writer(file.CreateWriter(windowBytes, counters));

    private sealed class Writer(MappedWindowWriter inner) : ISaSliceWriter
    {
        public void Write(long entryOffset, ReadOnlySpan<uint> data) => inner.WriteAt(entryOffset * 4, MemoryMarshal.AsBytes(data));

        public void Dispose() => inner.Dispose();
    }
}
