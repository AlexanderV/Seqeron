#pragma warning disable S2325, S1144, S1172, S3459, CS9113 // stub of the tests-first step: every member throws
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Ph5;

/// <summary>What a <see cref="SpillStack{T}"/> did with its spill file (ADR-004 section 3.2: "older blocks are appended to a spill file, read back backwards when a pop leaves the window").</summary>
internal readonly record struct SpillStats(long BlocksWritten, long BlocksRead, long BytesWritten, long BytesRead, long PeakFileBytes, long Flushes, long CleanEvictions);

/// <summary>
/// STUB (tests first): the external stack of the PH5 scan.
/// </summary>
internal sealed unsafe class SpillStack<T> : IDisposable
    where T : unmanaged
{
    internal SpillStack(BudgetedArena arena, string spillPath, int blockEntries, long flushBytes = 0) => throw new NotImplementedException();

    internal long Height => throw new NotImplementedException();

    internal long PeakHeight => throw new NotImplementedException();

    internal int BlockEntries => throw new NotImplementedException();

    internal ref T Top => throw new NotImplementedException();

    internal ref T TopForUpdate() => throw new NotImplementedException();

    internal void Push(in T value) => throw new NotImplementedException();

    internal void Pop() => throw new NotImplementedException();

    internal SpillStats Stats => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
