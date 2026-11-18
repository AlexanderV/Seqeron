#pragma warning disable S2325 // stub: the members become instance members of the implementation
namespace SuffixTree.Esa.External;

/// <summary>
/// The chunk arithmetic of PH4 (ADR-004 section 4.2): the text positions <c>0 .. n - 1</c> are cut into <see cref="TextChunks"/> chunks of R positions (the Phi tuples are distributed by
/// <c>SA[i] / R</c>), the ranks <c>0 .. n</c> into <see cref="RankChunks"/> chunks of R ranks (the LCP tuples are distributed by <c>i / R</c>). The division by R is the multiplication of <see cref="ChunkDivider"/>.
/// </summary>
internal readonly struct LcpGeometry
{
    internal LcpGeometry(long n, long chunkPositions) => throw new NotImplementedException();

    /// <summary>Text length n.</summary>
    internal long N => throw new NotImplementedException();

    /// <summary>R: positions per text chunk and ranks per rank chunk.</summary>
    internal long R => throw new NotImplementedException();

    /// <summary>Chunks of the text positions <c>0 .. n - 1</c> (zero for the empty text).</summary>
    internal long TextChunks => throw new NotImplementedException();

    /// <summary>Chunks of the ranks <c>0 .. n</c> (at least one).</summary>
    internal long RankChunks => throw new NotImplementedException();

    /// <summary>Positions in text chunk <paramref name="chunk"/>.</summary>
    internal long TextChunkLength(long chunk) => throw new NotImplementedException();

    /// <summary>Ranks in rank chunk <paramref name="chunk"/>.</summary>
    internal long RankChunkLength(long chunk) => throw new NotImplementedException();

    /// <summary><c>x / R</c> for a position or a rank.</summary>
    internal uint ChunkOf(uint x) => throw new NotImplementedException();
}
