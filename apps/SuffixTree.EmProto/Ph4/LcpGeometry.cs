using System.Runtime.CompilerServices;

namespace SuffixTree.EmProto.Ph4;

/// <summary>
/// The chunk arithmetic of PH4 (ADR-004 section 4.2): the text positions <c>0 .. n - 1</c> are cut into <see cref="TextChunks"/> chunks of R positions (the Phi tuples are distributed by
/// <c>SA[i] / R</c>), the ranks <c>0 .. n</c> into <see cref="RankChunks"/> chunks of R ranks (the LCP tuples are distributed by <c>i / R</c>). The division by R is done with a multiplication
/// (Lemire, Kaser, Kurz: "Faster remainder by direct computation", 2019: <c>floor(x / d) = high64(ceil(2^64 / d) * x)</c> for every 32-bit x and every d &gt;= 2), because it runs once per suffix-array entry.
/// </summary>
internal readonly struct LcpGeometry
{
    private readonly ulong _magic;

    internal LcpGeometry(long n, long chunkPositions)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        ArgumentOutOfRangeException.ThrowIfLessThan(chunkPositions, 1);
        if (n >= uint.MaxValue) throw new NotSupportedException("positions are 32-bit in the prototype");
        N = n;
        R = chunkPositions;
        TextChunks = (n + chunkPositions - 1) / chunkPositions;
        RankChunks = (n + chunkPositions) / chunkPositions; // ceil((n + 1) / R)
        _magic = chunkPositions is >= 2 and <= uint.MaxValue ? (ulong.MaxValue / (ulong)chunkPositions) + 1 : 0;
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
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal uint ChunkOf(uint x)
    {
        if (_magic != 0) return (uint)Math.BigMul(_magic, x, out _);
        return R == 1 ? x : 0; // R = 1, or R above 32 bits: a single chunk
    }
}
