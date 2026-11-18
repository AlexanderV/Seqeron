namespace SuffixTree.Esa.External;

/// <summary>
/// The chunk arithmetic of PH4 (ADR-004 section 4.2): the text positions <c>0 .. n - 1</c> are cut into <see cref="TextChunks"/> chunks of R positions (the Phi tuples are distributed by
/// <c>SA[i] / R</c>), the ranks <c>0 .. n</c> into <see cref="RankChunks"/> chunks of R ranks (the LCP tuples are distributed by <c>i / R</c>). The division by R is the multiplication of <see cref="ChunkDivider"/>,
/// because it runs once per suffix-array entry; R = 1 and R above 32 bits (a single chunk) are handled there.
/// </summary>
internal readonly struct LcpGeometry
{
    private readonly ChunkDivider _divide;

    /// <exception cref="NotSupportedException">n is above <see cref="PackedTextFormat.MaxSymbols"/> (positions are 32-bit).</exception>
    internal LcpGeometry(long n, long chunkPositions)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        ArgumentOutOfRangeException.ThrowIfLessThan(chunkPositions, 1);
        if (n > PackedTextFormat.MaxSymbols) throw new NotSupportedException($"{n} symbols: positions are 32-bit, the limit is {PackedTextFormat.MaxSymbols}");
        N = n;
        R = chunkPositions;
        TextChunks = (n + chunkPositions - 1) / chunkPositions;
        RankChunks = (n + chunkPositions) / chunkPositions; // ceil((n + 1) / R)
        _divide = new ChunkDivider(chunkPositions);
    }

    /// <summary>Text length n.</summary>
    internal long N { get; }

    /// <summary>R: positions per text chunk and ranks per rank chunk.</summary>
    internal long R { get; }

    /// <summary>Chunks of the text positions <c>0 .. n - 1</c> (zero for the empty text).</summary>
    internal long TextChunks { get; }

    /// <summary>Chunks of the ranks <c>0 .. n</c> (at least one).</summary>
    internal long RankChunks { get; }

    /// <summary>Positions in text chunk <paramref name="chunk"/>.</summary>
    internal long TextChunkLength(long chunk) => Math.Min(R, N - (chunk * R));

    /// <summary>Ranks in rank chunk <paramref name="chunk"/>.</summary>
    internal long RankChunkLength(long chunk) => Math.Min(R, N + 1 - (chunk * R));

    /// <summary><c>x / R</c> for a position or a rank.</summary>
    internal uint ChunkOf(uint x) => _divide.Of(x);
}
