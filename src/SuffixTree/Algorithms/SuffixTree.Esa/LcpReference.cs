namespace SuffixTree.Esa;

/// <summary>The 32-bit instantiation of <see cref="LcpReference{TIdx}"/> under the short name used by callers and tests.</summary>
internal static class LcpReference
{
    /// <summary>Builds the compact LCP array of the text from its suffix array (n + 1 layout, <c>sa[0] = n</c>).</summary>
    internal static CompactLcp Build<TReader>(TReader text, ReadOnlySpan<int> sa)
        where TReader : ISymbolReader, allows ref struct
        => LcpReference<int>.Build(text, sa);

    /// <summary>The permuted LCP array (length n), see <see cref="LcpReference{TIdx}.ComputePlcp"/>.</summary>
    internal static int[] ComputePlcp<TReader>(TReader text, ReadOnlySpan<int> sa)
        where TReader : ISymbolReader, allows ref struct
        => LcpReference<int>.ComputePlcp(text, sa);
}
