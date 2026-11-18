using System.Numerics;

namespace SuffixTree.Esa.External;

/// <summary>
/// The layout of the packed text, TEXT codec 4 (ADR-003 section 4 kind 1, ADR-004 section 4.2 row PH0, ADR-004 section 6a Q3): the dense symbols of the text (the rank of the code unit among the distinct code units, increasing)
/// stored with <c>b</c> bits per symbol, <c>b = max(1, ceil(log2 sigma))</c>, as one big-endian bit stream: symbol <c>i</c> occupies the stream bits <c>[i*b, (i+1)*b)</c>, the first bit of a byte is its most
/// significant bit, the unused bits of the last data byte are zero, and <see cref="PadBytes"/> zero bytes follow the data, so that a 64-bit window may start at any symbol up to the end (and the padding reads as
/// the symbol 0). This is the layout of the prototype (<c>apps/SuffixTree.EmProto</c>, <c>PackedText</c> and PH0 of E3-E8), which the pipeline of E6-E8 produced and <c>e8-verify</c> checked; the byte count of
/// the payload is <see cref="TotalBytes"/>. Nothing in the layout depends on <c>b</c> except the two arithmetic functions below: a wider code is a new codec id, not a change of this one.
/// </summary>
internal static class PackedTextFormat
{
    /// <summary>The codec id of the packed text in the TEXT section entry (ADR-003 reserves 4).</summary>
    internal const int CodecId = 4;

    /// <summary>Zero bytes after the last data byte (a 256-bit window plus the carry byte reads at most 41 bytes past the start byte).</summary>
    internal const int PadBytes = 64;

    /// <summary>Smallest bits per symbol (a text of one letter takes one bit).</summary>
    internal const int MinBits = 1;

    /// <summary>Largest bits per symbol of the file codec: byte alphabets (up to 256 distinct code units).</summary>
    internal const int MaxBits = 8;

    /// <summary>Design limit of the writer: the 32-bit positions of the builder (ADR-004 section 4.3, <c>w = 4</c>: n at most 2^32 - 4).</summary>
    internal const long MaxSymbols = uint.MaxValue - 3L;

    /// <summary>Smallest <c>b</c> with <c>2^b &gt;= sigma</c>, at least 1.</summary>
    internal static int BitsFor(long sigma)
    {
        if (sigma <= 2) return 1;
        return 64 - BitOperations.LeadingZeroCount((ulong)(sigma - 1));
    }

    /// <summary>Bytes of the data part for <paramref name="n"/> symbols of <paramref name="bits"/> bits (the last byte is padded with zero bits).</summary>
    internal static long DataBytes(long n, int bits) => ((n * bits) + 7) / 8;

    /// <summary>Bytes of the payload: the data part and the zero padding.</summary>
    internal static long TotalBytes(long n, int bits) => DataBytes(n, bits) + PadBytes;

    /// <summary>
    /// Packs dense symbols (each below <c>2^bits</c>) into <paramref name="destination"/>, which receives <c>DataBytes(symbols.Length, bits)</c> bytes (the unused bits of the last byte are written as zero; nothing beyond
    /// the data is touched). A block of a multiple of 8 symbols ends on a byte boundary, so consecutive blocks can be packed independently.
    /// </summary>
    internal static void PackInto(ReadOnlySpan<byte> symbols, int bits, Span<byte> destination)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(bits, MinBits);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(bits, MaxBits);
        long need = DataBytes(symbols.Length, bits);
        if (destination.Length < need) throw new ArgumentException($"the destination has {destination.Length} bytes but {symbols.Length} symbols of {bits} bits need {need}", nameof(destination));
        ulong acc = 0;
        int accBits = 0;
        int o = 0;
        foreach (byte s in symbols)
        {
            if ((s >> bits) != 0) throw new ArgumentOutOfRangeException(nameof(symbols), $"symbol {s} does not fit in {bits} bits");
            acc = (acc << bits) | s;
            accBits += bits;
            while (accBits >= 8)
            {
                destination[o++] = (byte)(acc >> (accBits - 8));
                accBits -= 8;
            }
        }
        if (accBits > 0) destination[o] = (byte)(acc << (8 - accBits));
    }
}

/// <summary>
/// What PH0 found about a text and what a reader needs to check one: the length, the width, the alphabet (the distinct code units in increasing order; the dense symbol <c>i</c> is <c>Letters[i]</c>), the symbol counts
/// <c>C[]</c> and the SHA-256 of the input bytes (lower-case hexadecimal), as the prototype's PH0 reported them. The fields map to the container as ADR-004 section 5a (1b.7b) proposes.
/// </summary>
internal sealed record PackedTextInfo(long N, int Bits, byte[] Letters, long[] Counts, string Sha256)
{
    /// <summary>Number of distinct code units.</summary>
    internal int Sigma => Letters.Length;

    /// <summary>Bytes of the data part.</summary>
    internal long DataBytes => PackedTextFormat.DataBytes(N, Bits);

    /// <summary>Bytes of the payload including the padding.</summary>
    internal long TotalBytes => PackedTextFormat.TotalBytes(N, Bits);
}
