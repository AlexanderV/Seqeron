#pragma warning disable S2325 // stub: the members become instance members of the implementation
namespace SuffixTree.Esa.External;

/// <summary>
/// What an <see cref="AlignedWindowWriter"/> did: windows mapped (<see cref="Opened"/>), windows mapped again after having been released (<see cref="Reopened"/>: a page that was flushed and is written again is written to storage
/// a second time, the write amplification of layout A), flushes (msync) and the bytes of the windows they covered, 32-bit writes, windows released to make room, distinct windows ever mapped.
/// </summary>
internal readonly record struct WindowStats(long Opened, long Reopened, long Flushes, long FlushedBytes, long Writes, long Evictions, long DistinctWindows);

/// <summary>
/// Writes 32-bit values at arbitrary entries of a range of a <see cref="MappedOutputFile"/> through a small set of aligned mapped windows (the child-table fields of layout A are written behind the cursor "through bounded
/// windows", ADR-004 section 3.2): a write maps the window <c>offset / W</c> (aligned, at most <c>slots</c> windows at a time, the least recently used released first); with <c>flushOnRelease</c> a released window is
/// flushed (msync), so the dirty pages of the file stay within <c>slots x W</c> plus what the kernel has not written yet. The file is created by the caller (<see cref="MappedOutputFile.Create"/>) and is not owned
/// here. Not thread-safe.
/// </summary>
internal sealed class AlignedWindowWriter : IDisposable
{
    /// <param name="file">The output file; this writer counts itself as one of its writers until it is disposed.</param>
    /// <param name="start">Offset of the range in the file.</param>
    /// <param name="length">Length of the range in bytes (a multiple of 4).</param>
    /// <param name="windowBytes">Size of a window: a power of two of at least 4.</param>
    /// <param name="slots">Windows mapped at a time.</param>
    /// <param name="flushOnRelease">msync a window when it is released.</param>
    /// <param name="counters">Receives the windows and the flushed bytes.</param>
    internal AlignedWindowWriter(MappedOutputFile file, long start, long length, int windowBytes, int slots, bool flushOnRelease, IoCounters? counters = null) => throw new NotImplementedException();

    /// <summary>Length of the range in bytes.</summary>
    internal long Length => throw new NotImplementedException();

    internal WindowStats Stats => throw new NotImplementedException();

    /// <summary>Writes the 32-bit little-endian value at entry <paramref name="index"/> (byte offset <c>4 x index</c> of the range).</summary>
    internal void Write32(long index, uint value) => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
