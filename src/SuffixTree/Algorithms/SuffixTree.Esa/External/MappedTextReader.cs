#pragma warning disable // STUB (tests first): removed by the implementation commit
namespace SuffixTree.Esa.External;

/// <summary>STUB (tests first): long-indexed random reads of a b-bit packed text (TEXT codec 4) through a memory-mapped file.</summary>
internal sealed unsafe class MappedTextReader : IDisposable
{
    private MappedTextReader()
    {
    }

    internal long Length => throw new NotImplementedException();

    internal int Bits => throw new NotImplementedException();

    internal int KeySymbols => throw new NotImplementedException();

    internal int KeyShift => throw new NotImplementedException();

    internal int VectorSymbols => throw new NotImplementedException();

    internal bool UseVectors
    {
        get => throw new NotImplementedException();
        set => throw new NotImplementedException();
    }

    internal static MappedTextReader Open(string path, long n, int bits) => throw new NotImplementedException();

    internal static MappedTextReader OpenRegion(string path, long offset, long length, long n, int bits) => throw new NotImplementedException();

    internal static MappedTextReader FromMemory(byte* data, long length, long n, int bits) => throw new NotImplementedException();

    internal ulong Window(long i) => throw new NotImplementedException();

    internal ulong Key(long i) => throw new NotImplementedException();

    internal uint Symbol(long i) => throw new NotImplementedException();

    internal void ReadSymbols(long start, Span<byte> destination) => throw new NotImplementedException();

    internal long Lce(long i, long j, long max, ref long steps) => throw new NotImplementedException();

    internal int CompareSuffixes(long i, long j, long skip, ref long steps, out long lce) => throw new NotImplementedException();

    internal void Verify(PackedTextInfo info, bool checkHash = true) => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
