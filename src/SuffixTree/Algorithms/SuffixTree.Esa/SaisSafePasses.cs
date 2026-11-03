namespace SuffixTree.Esa;

/// <summary>
/// The induction passes and the naming of <see cref="SaisKernel"/>, verbatim copies of the algorithm of <see cref="SaisReference{TIdx}"/> in safe code (Span
/// indexers with bounds checks), specialised to <c>int</c> and to the byte-per-position type array of the kernel. They are the part of SA-IS that iteration 1o.1b
/// deliberately leaves alone: iteration 1o.1c replaces them by pointer loops. They exist as a copy (and not as calls into the reference) so that the permanent oracle
/// stays literally untouched. The bucket array is prepared by the caller (<c>SetBucketHeads</c> before <see cref="InduceL{TReader}"/>, <c>SetBucketTails</c> before <see cref="InduceS{TReader}"/>).
/// </summary>
internal static class SaisSafePasses
{
    /// <summary>Induces the L-type suffixes from the sorted entries, left to right.</summary>
    internal static void InduceL<TReader>(TReader t, Span<int> sa, ReadOnlySpan<byte> types, Span<int> bkt)
        where TReader : ISymbolReader, allows ref struct
    {
        for (int i = 0; i <= t.Length; i++)
        {
            int p = sa[i];
            if (p <= 0) continue; // empty slot, or suffix 0 which has no predecessor
            int j = p - 1;
            if (types[j] != 0) continue;
            int c = t[j];
            int slot = bkt[c];
            sa[slot] = j;
            bkt[c] = slot + 1;
        }
    }

    /// <summary>Induces the S-type suffixes from the sorted entries, right to left.</summary>
    internal static void InduceS<TReader>(TReader t, Span<int> sa, ReadOnlySpan<byte> types, Span<int> bkt)
        where TReader : ISymbolReader, allows ref struct
    {
        for (int i = t.Length; i >= 1; i--)
        {
            int p = sa[i];
            if (p <= 0) continue;
            int j = p - 1;
            if (types[j] == 0) continue;
            int c = t[j];
            int slot = bkt[c] - 1;
            sa[slot] = j;
            bkt[c] = slot;
        }
    }

    /// <summary>Names the sorted LMS substrings (equal substrings share a name) and stores the name of the LMS position p at <c>SA[1 + m + (p &gt;&gt; 1)]</c>. Returns the number of distinct names.</summary>
    internal static int NameLmsSubstrings<TReader>(TReader t, Span<int> sa, ReadOnlySpan<byte> types, int m)
        where TReader : ISymbolReader, allows ref struct
    {
        sa.Slice(m + 1).Fill(-1);
        int name = -1;
        int previous = -1;
        for (int i = 0; i < m; i++)
        {
            int pos = sa[1 + i];
            if (previous < 0 || !LmsSubstringsEqual(t, types, previous, pos)) name++;
            sa[1 + m + (pos >> 1)] = name;
            previous = pos;
        }
        return name + 1;
    }

    /// <summary>Equality of the LMS substrings that start at the distinct LMS positions a and b (same length, same symbols, same types); the one that reaches the sentinel is unique.</summary>
    internal static bool LmsSubstringsEqual<TReader>(TReader t, ReadOnlySpan<byte> types, int a, int b)
        where TReader : ISymbolReader, allows ref struct
    {
        int n = t.Length;
        for (int d = 0; a + d <= n && b + d <= n; d++)
        {
            int x = a + d, y = b + d;
            if (x == n || y == n) return false;
            if (t[x] != t[y] || types[x] != types[y]) return false;
            if (d > 0 && x > 0 && types[x] != 0 && types[x - 1] == 0) return true;
        }
        return false;
    }
}
