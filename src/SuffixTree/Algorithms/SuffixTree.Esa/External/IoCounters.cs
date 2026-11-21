namespace SuffixTree.Esa.External;

/// <summary>Values of an <see cref="IoCounters"/> at one moment.</summary>
internal readonly record struct IoSnapshot(
    long ReadBytes,
    long WriteBytes,
    long ReadCalls,
    long WriteCalls,
    long MappedWriteBytes,
    long MappedWindows,
    long MappedMaxWindowBytes,
    long ScratchLiveBytes,
    long ScratchPeakBytes)
{
    /// <summary>Explicit reads and writes plus the writes through mapped windows.</summary>
    internal long TotalBytes => ReadBytes + WriteBytes + MappedWriteBytes;

    /// <summary>Difference of the cumulative fields; the maxima and the live scratch are those of this snapshot (a maximum or a level is not a difference).</summary>
    internal IoSnapshot Since(IoSnapshot earlier)
        => new(
            ReadBytes - earlier.ReadBytes,
            WriteBytes - earlier.WriteBytes,
            ReadCalls - earlier.ReadCalls,
            WriteCalls - earlier.WriteCalls,
            MappedWriteBytes - earlier.MappedWriteBytes,
            MappedWindows - earlier.MappedWindows,
            MappedMaxWindowBytes,
            ScratchLiveBytes,
            ScratchPeakBytes);
}

/// <summary>
/// Byte and call counters of the explicit file I/O (every <see cref="CountedFile"/> that was given the instance), of the writes through mapped windows and of the scratch files of one build (live and peak bytes: ADR-004 Q1
/// reports peak scratch next to the total I/O). Thread-safe; one instance per build, handed to the pieces through their options (there is no process-wide state). Accesses through memory mappings that are not writes of
/// a <see cref="MappedWindowWriter"/> are not counted: the block-layer view of those is <see cref="ProcIoSnapshot"/>.
/// </summary>
internal sealed class IoCounters
{
    private long _readBytes;
    private long _writeBytes;
    private long _readCalls;
    private long _writeCalls;
    private long _mappedWriteBytes;
    private long _mappedWindows;
    private long _mappedMaxWindow;
    private long _scratchLive;
    private long _scratchPeak;
    private long _flushes;
    private long _reserved;

    internal long ReadBytes => Interlocked.Read(ref _readBytes);

    internal long WriteBytes => Interlocked.Read(ref _writeBytes);

    internal long ReadCalls => Interlocked.Read(ref _readCalls);

    internal long WriteCalls => Interlocked.Read(ref _writeCalls);

    internal long MappedWriteBytes => Interlocked.Read(ref _mappedWriteBytes);

    /// <summary>Windows (views) created.</summary>
    internal long MappedWindows => Interlocked.Read(ref _mappedWindows);

    /// <summary>Largest window created, in bytes.</summary>
    internal long MappedMaxWindowBytes => Interlocked.Read(ref _mappedMaxWindow);

    /// <summary>Bytes of scratch files alive now.</summary>
    internal long ScratchLiveBytes => Interlocked.Read(ref _scratchLive);

    /// <summary>The largest <see cref="ScratchLiveBytes"/> that was ever reached.</summary>
    internal long ScratchPeakBytes => Interlocked.Read(ref _scratchPeak);

    /// <summary>Forced writes to storage (fsync) through <see cref="CountedFile.Flush"/>: the calls, not the bytes (an fsync has no byte count; the order of the two fsyncs of the publication of a container is observable here).</summary>
    internal long Flushes => Interlocked.Read(ref _flushes);

    internal void AddRead(long bytes)
    {
        Interlocked.Add(ref _readBytes, bytes);
        Interlocked.Increment(ref _readCalls);
    }

    internal void AddWrite(long bytes)
    {
        Interlocked.Add(ref _writeBytes, bytes);
        Interlocked.Increment(ref _writeCalls);
    }

    internal void AddFlush() => Interlocked.Increment(ref _flushes);

    /// <summary>Bytes of zeros written to reserve the blocks of a range that was added to an existing file (<see cref="MappedOutputFile.Grow"/>): reported on their own, not part of <see cref="WriteBytes"/>, <see cref="IoSnapshot.TotalBytes"/> or the snapshot (the process's <c>wchar</c> has them).</summary>
    internal long ReservedWriteBytes => Interlocked.Read(ref _reserved);

    internal void AddReserved(long bytes) => Interlocked.Add(ref _reserved, bytes);

    internal void AddMappedWrite(long bytes) => Interlocked.Add(ref _mappedWriteBytes, bytes);

    internal void AddWindow(long bytes)
    {
        Interlocked.Increment(ref _mappedWindows);
        Atomic.RaiseTo(ref _mappedMaxWindow, bytes);
    }

    internal void AddScratch(long bytes) => Atomic.RaiseTo(ref _scratchPeak, Interlocked.Add(ref _scratchLive, bytes));

    internal void RemoveScratch(long bytes) => Interlocked.Add(ref _scratchLive, -bytes);

    internal void Reset()
    {
        Interlocked.Exchange(ref _readBytes, 0);
        Interlocked.Exchange(ref _writeBytes, 0);
        Interlocked.Exchange(ref _readCalls, 0);
        Interlocked.Exchange(ref _writeCalls, 0);
        Interlocked.Exchange(ref _mappedWriteBytes, 0);
        Interlocked.Exchange(ref _mappedWindows, 0);
        Interlocked.Exchange(ref _mappedMaxWindow, 0);
        Interlocked.Exchange(ref _scratchLive, 0);
        Interlocked.Exchange(ref _scratchPeak, 0);
        Interlocked.Exchange(ref _flushes, 0);
        Interlocked.Exchange(ref _reserved, 0);
    }

    /// <summary>The values now (each read atomically, the set as a whole is not a transaction).</summary>
    internal IoSnapshot Snapshot()
        => new(ReadBytes, WriteBytes, ReadCalls, WriteCalls, MappedWriteBytes, MappedWindows, MappedMaxWindowBytes, ScratchLiveBytes, ScratchPeakBytes);
}
