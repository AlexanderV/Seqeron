namespace SuffixTree.EmProto.Gen;

/// <summary>What a streamed generation produced: the length, the seed, the SHA-256 of the whole output and the count of every letter of <see cref="GenomeStream.Letters"/>.</summary>
internal sealed record GenomeStreamResult(long N, ulong Seed, string Sha256, long[] Counts, double Seconds);

/// <summary>
/// E8 (ADR-004 section 5, row E8): the genome-like text of <see cref="Families.GenomeLike"/> produced as a stream of ASCII bytes with a bounded buffer, for any n (1 to 3x10^9 and beyond), byte for byte equal to
/// <c>Families.ToAscii(Families.GenomeLike(n, seed))</c> wherever that can be built in RAM, and identical for any buffer size.
/// </summary>
internal static class GenomeStream
{
    /// <summary>The seed of the family in <c>FamilyLoader</c> (<c>genome:N</c>), i.e. of the inputs of E3 to E7.</summary>
    internal const ulong DefaultSeed = 20240002;

    /// <summary>The letters of the codes 0..4 (A C G N T), as <see cref="Families.GenomeLike"/>.</summary>
    internal const string Letters = "ACGNT";

    /// <summary>Streams the text to <paramref name="sink"/> in blocks of at most <paramref name="bufferBytes"/> bytes.</summary>
    internal static GenomeStreamResult Generate(long n, ulong seed, Action<ReadOnlySpan<byte>> sink, int bufferBytes = 1 << 20) => throw new NotImplementedException();

    /// <summary>The text of the checksum file (<c>key=value</c> lines).</summary>
    internal static string SumsText(GenomeStreamResult r) => throw new NotImplementedException();
}
