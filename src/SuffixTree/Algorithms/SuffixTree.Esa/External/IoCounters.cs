#pragma warning disable S2325, CA1822, S1144, IDE0051, IDE0060, S1172 // red stub of work package 1b.7a: every member throws NotImplementedException; the implementation replaces this file

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
    internal long TotalBytes => throw new NotImplementedException();

    /// <summary>Difference of the cumulative fields; the maxima and the live scratch are those of this snapshot.</summary>
    internal IoSnapshot Since(IoSnapshot earlier) => throw new NotImplementedException();
}

/// <summary>Byte and call counters of the explicit file I/O, the writes through mapped windows and the scratch files of one build. Thread-safe.</summary>
internal sealed class IoCounters
{
    internal long ReadBytes => throw new NotImplementedException();

    internal long WriteBytes => throw new NotImplementedException();

    internal long ReadCalls => throw new NotImplementedException();

    internal long WriteCalls => throw new NotImplementedException();

    internal long MappedWriteBytes => throw new NotImplementedException();

    internal long MappedWindows => throw new NotImplementedException();

    internal long MappedMaxWindowBytes => throw new NotImplementedException();

    internal long ScratchLiveBytes => throw new NotImplementedException();

    internal long ScratchPeakBytes => throw new NotImplementedException();

    internal void AddRead(long bytes) => throw new NotImplementedException();

    internal void AddWrite(long bytes) => throw new NotImplementedException();

    internal void AddMappedWrite(long bytes) => throw new NotImplementedException();

    internal void AddWindow(long bytes) => throw new NotImplementedException();

    internal void AddScratch(long bytes) => throw new NotImplementedException();

    internal void RemoveScratch(long bytes) => throw new NotImplementedException();

    internal void Reset() => throw new NotImplementedException();

    internal IoSnapshot Snapshot() => throw new NotImplementedException();
}
