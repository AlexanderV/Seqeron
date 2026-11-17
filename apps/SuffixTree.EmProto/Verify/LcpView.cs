namespace SuffixTree.EmProto.Verify;

/// <summary>
/// Read-only view of the compact LCP of an <c>e6</c> run: the base array of 1, 2 or 4 bytes per rank (little endian) and the overflow pairs <c>(rank, value)</c> of two 4-byte fields, the pairs in ascending rank order, one per base
/// entry that holds the escape value (255, 65,535; the 4-byte base has none). The pointers belong to the caller (mapped files).
/// </summary>
internal sealed unsafe class LcpView(byte* basePointer, int width, uint* pairs, long pairCount)
{
    internal int Width { get; } = width;

    /// <summary>The value that marks "look in the overflow list" (none for the 4-byte base).</summary>
    internal uint Escape { get; } = width switch { 1 => byte.MaxValue, 2 => ushort.MaxValue, _ => uint.MaxValue };

    internal long PairCount { get; } = pairCount;

    internal bool HasEscape => Width < 4;

    internal uint RankOfPair(long k) => pairs[2 * k];

    internal uint ValueOfPair(long k) => pairs[(2 * k) + 1];

    /// <summary>The stored base entry (the escape value for an overflow entry).</summary>
    internal uint BaseAt(long rank) => Width switch
    {
        1 => basePointer[rank],
        2 => ((ushort*)basePointer)[rank],
        _ => ((uint*)basePointer)[rank],
    };

    /// <summary>Index of the first pair whose rank is at least <paramref name="rank"/> (the pair count when there is none).</summary>
    internal long LowerBound(long rank)
    {
        long lo = 0, hi = PairCount;
        while (lo < hi)
        {
            long mid = lo + ((hi - lo) >> 1);
            if (pairs[2 * mid] < rank) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }
}
