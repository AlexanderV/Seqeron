namespace SuffixTree.Esa;

/// <summary>
/// Hypothesis H-OCC (ADR-001 section 4): the edge of every internal node of the Ukkonen tree records the LEFTMOST occurrence of the node's string,
/// <c>min SA[lb..rb]</c>. <c>Traverse</c> needs it to reproduce the edge positions of the Ukkonen tree (they are part of
/// <c>SuffixTreeSerializer.CalculateLogicalHash</c>). One bottom-up scan over the LCP stream computes it for every internal interval.
/// </summary>
internal static class LeftmostOccurrence
{
    /// <summary>
    /// The leftmost occurrence of every internal interval, stored at the rank of the interval's FIRST l-index (that rank belongs to exactly one interval);
    /// all other entries are 0. The array has <c>n + 1</c> entries. Transient: 4 bytes per rank, used only by the diagnostic traversal.
    /// </summary>
    internal static int[] Compute(int[] sa, CompactLcp lcp)
    {
        ArgumentNullException.ThrowIfNull(sa);
        ArgumentNullException.ThrowIfNull(lcp);
        int n = sa.Length - 1;
        if (lcp.Length != n + 1) throw new ArgumentException("the LCP array and the suffix array differ in length", nameof(lcp));
        var occurrence = new int[n + 1];
        if (n == 0) return occurrence;

        int capacity = 64;
        int[] ell = new int[capacity], first = new int[capacity], min = new int[capacity];
        int sp = 1;
        first[0] = -1;
        min[0] = int.MaxValue;
        var stream = new CompactLcpStream(lcp);
        for (int i = 1; i <= n + 1; i++)
        {
            int cur = i <= n ? stream.Next() : -1;
            // the leaf of rank i - 1 is complete: it belongs to the innermost open interval (a new interval that starts at it is pushed below and inherits it)
            if (sa[i - 1] < min[sp - 1]) min[sp - 1] = sa[i - 1];
            int carry = sa[i - 1];
            while (sp > 0 && ell[sp - 1] > cur)
            {
                sp--;
                occurrence[first[sp]] = min[sp];
                carry = min[sp];
                if (sp > 0 && carry < min[sp - 1]) min[sp - 1] = carry;
            }
            if (i > n) break;

            int top = sp - 1;
            if (ell[top] < cur)
            {
                if (sp == capacity)
                {
                    capacity *= 2;
                    Array.Resize(ref ell, capacity);
                    Array.Resize(ref first, capacity);
                    Array.Resize(ref min, capacity);
                }
                ell[sp] = cur;
                first[sp] = i;
                min[sp] = carry;
                sp++;
            }
            else if (first[top] < 0)
            {
                first[top] = i;
            }
        }
        return occurrence;
    }
}
