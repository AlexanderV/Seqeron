#pragma warning disable // STUB (tests first): removed by the implementation commit
namespace SuffixTree.Esa.External;

/// <summary>STUB (tests first): layout constants and arithmetic of the packed text, TEXT codec 4 (ADR-003 kind 1, ADR-004 PH0).</summary>
internal static class PackedTextFormat
{
    /// <summary>The codec id of the packed text in the TEXT section entry.</summary>
    internal const int CodecId = 4;

    /// <summary>Zero bytes after the last data byte.</summary>
    internal const int PadBytes = 64;

    /// <summary>Smallest bits per symbol.</summary>
    internal const int MinBits = 1;

    /// <summary>Largest bits per symbol of the file codec (byte alphabets).</summary>
    internal const int MaxBits = 8;

    /// <summary>Design limit of the writer: the 32-bit positions of the builder.</summary>
    internal const long MaxSymbols = uint.MaxValue - 3L;

    internal static int BitsFor(long sigma) => throw new NotImplementedException();

    internal static long DataBytes(long n, int bits) => throw new NotImplementedException();

    internal static long TotalBytes(long n, int bits) => throw new NotImplementedException();

    internal static void PackInto(ReadOnlySpan<byte> symbols, int bits, Span<byte> destination) => throw new NotImplementedException();
}

/// <summary>What PH0 found about a text and what a reader needs to check one: length, width, alphabet, symbol counts, SHA-256 of the input bytes (lower-case hex).</summary>
internal sealed record PackedTextInfo(long N, int Bits, byte[] Letters, long[] Counts, string Sha256)
{
    internal int Sigma => Letters.Length;

    internal long DataBytes => PackedTextFormat.DataBytes(N, Bits);

    internal long TotalBytes => PackedTextFormat.TotalBytes(N, Bits);
}
