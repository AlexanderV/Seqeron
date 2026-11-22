// STUB (tests first): the signatures of the hand-off query set, the streamed oracle and the suffix-array search. Replaced by the implementation in the next commit.
namespace SuffixTree.Esa.Bench.Links;

internal sealed record HandoffQuery(int Index, bool Present, long Origin, byte[] Bytes);

internal sealed record ScanResult(long[] Counts, List<long>[] Positions);

internal static class HandoffQuerySet
{
    internal static byte[] AlphabetOf(ReadOnlySpan<byte> sample) => throw new NotImplementedException();

    internal static IReadOnlyList<HandoffQuery> Build(long n, Func<long, int, byte[]> read, ReadOnlySpan<byte> alphabet, IReadOnlyList<int> lengths, int present, int absent, ulong seed) => throw new NotImplementedException();
}

internal static class RollingScan
{
    internal static ScanResult Run(Stream text, IReadOnlyList<HandoffQuery> queries, int cap, int blockBytes = 4 << 20) => throw new NotImplementedException();
}

internal static class SaSearch
{
    internal static (long Lo, long Hi) Range(ReadOnlySpan<byte> query, long n, Func<long, long> sa, Func<long, byte[], int> readText) => throw new NotImplementedException();
}
