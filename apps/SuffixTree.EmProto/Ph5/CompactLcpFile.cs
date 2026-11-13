#pragma warning disable S2325, S1144, S1172, S3459, CS9113 // stub of the tests-first step: every member throws
namespace SuffixTree.EmProto.Ph5;

/// <summary>Receives blocks of resolved LCP values in the order a file reader produces them.</summary>
internal delegate void LcpBlockHandler(long firstRank, ReadOnlySpan<uint> values);

/// <summary>STUB (tests first): the compact LCP of PH4 (base file of 1, 2 or 4 bytes per entry plus overflow pairs) read in rank order or backwards.</summary>
internal sealed unsafe class CompactLcpFile : IDisposable
{
    internal CompactLcpFile(string basePath, string overflowPath, long n, int bufferBytes = 1 << 20) => throw new NotImplementedException();

    internal int Width => throw new NotImplementedException();

    internal long OverflowCount => throw new NotImplementedException();

    internal void ForEachForward(LcpBlockHandler handler) => throw new NotImplementedException();

    internal void ForEachBackward(LcpBlockHandler handler) => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
