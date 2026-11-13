#pragma warning disable S2325, S1144, S1172, S3459, CS9113 // stub of the tests-first step: every member throws
namespace SuffixTree.EmProto.Ph5;

/// <summary>What an <see cref="AlignedWindowWriter"/> did: windows mapped, windows mapped again after having been released, flushes (msync) and the bytes they covered, writes.</summary>
internal readonly record struct WindowStats(long Opened, long Reopened, long Flushes, long FlushedBytes, long Writes, long Evictions, long DistinctWindows);

/// <summary>STUB (tests first): a file of known length written through a small set of aligned mapped windows.</summary>
internal sealed unsafe class AlignedWindowWriter : IDisposable
{
    internal AlignedWindowWriter(string path, long length, int windowBytes, int slots, bool flushOnRelease) => throw new NotImplementedException();

    internal long Length => throw new NotImplementedException();

    internal void Write32(long index, uint value) => throw new NotImplementedException();

    internal void Write(long offset, ReadOnlySpan<byte> data) => throw new NotImplementedException();

    internal void Flush() => throw new NotImplementedException();

    internal WindowStats Stats => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
