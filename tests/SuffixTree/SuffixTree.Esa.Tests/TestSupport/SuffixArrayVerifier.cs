namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>
/// Linear-time validity check of a suffix array (n + 1 layout), independent of any construction algorithm and of the naive
/// sort: <c>sa</c> is a permutation of 0..n with <c>sa[0] = n</c>, and for every adjacent pair <c>i = sa[r-1]</c>, <c>j = sa[r]</c>
/// either <c>text[i] &lt; text[j]</c> or (<c>text[i] == text[j]</c> and <c>rank[i+1] &lt; rank[j+1]</c>), where <c>rank[n] = 0</c>.
/// By induction on the suffix length this is equivalent to "sa is the sorted order of the suffixes" (ordinal code-unit order).
/// Used where the naive O(n^2 log n) oracle is too slow (n up to 10^6); cross-validated against the naive oracle on small texts.
/// </summary>
public static class SuffixArrayVerifier
{
    /// <summary>Returns null when <paramref name="sa"/> is the suffix array of <paramref name="text"/>, otherwise a description of the first violation.</summary>
    public static string? Check(string text, ReadOnlySpan<int> sa)
    {
        int n = text.Length;
        if (sa.Length != n + 1) return $"length {sa.Length} != n + 1 = {n + 1}";
        if (sa[0] != n) return $"sa[0] = {sa[0]} != n = {n}";
        var rank = new int[n + 1];
        Array.Fill(rank, -1);
        for (int r = 0; r <= n; r++)
        {
            int p = sa[r];
            if ((uint)p > (uint)n) return $"sa[{r}] = {p} out of range";
            if (rank[p] != -1) return $"position {p} appears twice (ranks {rank[p]} and {r})";
            rank[p] = r;
        }
        for (int r = 2; r <= n; r++)
        {
            int i = sa[r - 1], j = sa[r];
            // i, j < n here because rank 0 is the empty suffix and the positions are distinct non-empty suffixes
            if (text[i] < text[j]) continue;
            if (text[i] > text[j]) return $"rank {r}: text[{i}] > text[{j}]";
            if (rank[i + 1] >= rank[j + 1]) return $"rank {r}: suffixes {i} and {j} start equal but their successors are out of order";
        }
        return null;
    }
}
