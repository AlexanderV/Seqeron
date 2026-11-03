namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>
/// Definitional suffix array (test oracle L0). Layout of ADR-002 section 2: <c>n + 1</c> entries, <c>SA[0] = n</c> (the empty
/// suffix, i.e. the virtual terminator that is smaller than every symbol), <c>SA[1..n]</c> the non-empty suffixes in
/// ascending ordinal UTF-16 code-unit order (a proper prefix sorts before its extensions).
/// </summary>
public static class NaiveSuffixArray
{
    /// <summary>Builds the suffix array by sorting the suffix start positions with ordinal span comparison.</summary>
    public static int[] Build(string text)
    {
        ArgumentNullException.ThrowIfNull(text);
        int n = text.Length;
        var sa = new int[n + 1];
        sa[0] = n;
        for (int i = 0; i < n; i++) sa[i + 1] = i;
        Array.Sort(sa, 1, n, Comparer<int>.Create((a, b) => text.AsSpan(a).SequenceCompareTo(text.AsSpan(b))));
        return sa;
    }
}

/// <summary>
/// Definitional LCP array (test oracle L0): <c>LCP[r]</c> is the length of the longest common prefix of the suffixes
/// <c>SA[r-1]</c> and <c>SA[r]</c> for <c>r &gt;= 1</c>; <c>LCP[0] = 0</c> and <c>LCP[1] = 0</c> (the predecessor of rank 1 is the empty suffix).
/// Length <c>n + 1</c>.
/// </summary>
public static class NaiveLcp
{
    /// <summary>Builds the LCP array by comparing every pair of adjacent suffixes directly.</summary>
    public static int[] Build(string text, int[] sa)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(sa);
        if (sa.Length != text.Length + 1) throw new ArgumentException("sa must have text.Length + 1 entries", nameof(sa));
        var lcp = new int[sa.Length];
        for (int r = 2; r < sa.Length; r++)
            lcp[r] = text.AsSpan(sa[r - 1]).CommonPrefixLength(text.AsSpan(sa[r]));
        return lcp;
    }
}

/// <summary>Kasai, Lee, Arimura, Arikawa, Park (2001) linear-time LCP from the inverse suffix array; same layout as <see cref="NaiveLcp"/>.</summary>
public static class KasaiLcp
{
    /// <summary>Builds the LCP array in O(n) from the text and its suffix array.</summary>
    public static int[] Build(string text, int[] sa)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(sa);
        int n = text.Length;
        if (sa.Length != n + 1) throw new ArgumentException("sa must have text.Length + 1 entries", nameof(sa));
        var rank = new int[n + 1];
        for (int r = 0; r <= n; r++) rank[sa[r]] = r;
        var lcp = new int[n + 1];
        int h = 0;
        for (int i = 0; i < n; i++)
        {
            int r = rank[i];
            if (r == 1)
            {
                h = 0; // predecessor is the empty suffix
                continue;
            }
            int j = sa[r - 1];
            while (i + h < n && j + h < n && text[i + h] == text[j + h]) h++;
            lcp[r] = h;
            if (h > 0) h--;
        }
        return lcp;
    }
}
