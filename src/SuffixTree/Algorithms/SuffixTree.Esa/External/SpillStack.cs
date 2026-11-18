#pragma warning disable S2325, CA1822, S1144, IDE0051, IDE0060, S1172 // red stub of work package 1b.7a: every member throws NotImplementedException; the implementation replaces this file

namespace SuffixTree.Esa.External;

/// <summary>What a <see cref="SpillStack{T}"/> did with its spill file.</summary>
internal readonly record struct SpillStats(long BlocksWritten, long BlocksRead, long BytesWritten, long BytesRead, long PeakFileBytes, long Flushes, long CleanEvictions);

/// <summary>The external stack of PH5 (ADR-004 section 3.2): two block buffers in RAM, older blocks spilled to a file, read back when a pop leaves the window. Not thread-safe.</summary>
internal sealed unsafe class SpillStack<T> : IDisposable
    where T : unmanaged
{
    internal SpillStack(BudgetedArena arena, string spillPath, int blockEntries, long flushBytes = 0, IoCounters? counters = null) => throw new NotImplementedException();

    internal long Height => throw new NotImplementedException();

    internal long PeakHeight => throw new NotImplementedException();

    internal int BlockEntries => throw new NotImplementedException();

    internal SpillStats Stats => throw new NotImplementedException();

    internal ref T Top => throw new NotImplementedException();

    internal ref T TopForUpdate() => throw new NotImplementedException();

    internal void Push(in T value) => throw new NotImplementedException();

    internal void Pop() => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
