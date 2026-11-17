namespace SuffixTree.EmProto.Verify;

/// <summary>
/// E8 of ADR-004 section 5: the checker of the outputs of <c>e6</c> (packed text, suffix array, compact LCP, child table, statistics) that needs no in-RAM reference and so works above 2^31 symbols. See <see cref="Names"/> for the checks.
/// </summary>
internal static class OutputVerifier
{
    /// <summary>The checks in the order they run.</summary>
    internal static readonly string[] Names = ["files", "meta", "text", "input", "sa_permutation", "sa_order", "lcp_values", "lcp_width", "stats", "child", "sha256"];

    internal static VerifyResult Run(VerifyOptions options, Action<string>? progress = null) => throw new NotImplementedException();
}
