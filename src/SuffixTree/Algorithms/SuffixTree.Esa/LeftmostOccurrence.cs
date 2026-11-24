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
        return Compute(new SaView(sa), lcp);
    }

    /// <summary>
    /// The same for a suffix array that is a view (an index that reads a container), at 64-bit ranks, positions and depths over any stored width: the scan reads the views through their long members. The table is a managed array
    /// of 4-byte positions, so a text whose <c>n + 1</c> ranks do not fit an array has no such table.
    /// </summary>
    /// <exception cref="NotSupportedException">The text has more than <see cref="Array.MaxLength"/> - 1 symbols: its table of occurrences cannot exist (the traversal of such a text is not supported).</exception>
    /// <exception cref="InvalidDataException">The LCP array is inconsistent (a negative value, or a value that leaves an interval without its first l-index).</exception>
    internal static int[] Compute(SaView sa, CompactLcp lcp)
    {
        ArgumentNullException.ThrowIfNull(lcp);
        long n = sa.Length - 1;
        if (n + 1 > Array.MaxLength) throw new NotSupportedException($"the table of the leftmost occurrences of a text of {n} symbols has {n + 1} entries, more than an array holds ({Array.MaxLength}): the traversal of a text this long is not supported (LeftmostOccurrence.Compute)");
        if (lcp.Length != n + 1) throw new ArgumentException("the LCP array and the suffix array differ in length", nameof(lcp));
        var occurrence = new int[n + 1];
        if (n == 0) return occurrence;

        int capacity = 64;
        long[] ell = new long[capacity], first = new long[capacity], min = new long[capacity];
        int sp = 1;
        first[0] = -1;
        min[0] = long.MaxValue;
        CompactLcp.LongCursor stream = lcp.CreateLongCursor();
        for (long i = 1; i <= n + 1; i++)
        {
            long cur = i <= n ? stream.Get(i) : -1;
            // the leaf of rank i - 1 is complete: it belongs to the innermost open interval (a new interval that starts at it is pushed below and inherits it)
            long leaf = sa[i - 1];
            if (leaf < min[sp - 1]) min[sp - 1] = leaf;
            long carry = leaf;
            while (sp > 0 && ell[sp - 1] > cur)
            {
                sp--;
                if (first[sp] < 0) throw DamagedLcp(); // only an LCP[1] other than 0 leaves the base interval without a first l-index
                occurrence[first[sp]] = IntBridge.Value(min[sp], "a leftmost occurrence (a position of a text below 2^31 symbols)");
                carry = min[sp];
                if (sp > 0 && carry < min[sp - 1]) min[sp - 1] = carry;
            }
            if (i > n) break;
            if (sp == 0) throw DamagedLcp(); // a negative LCP value popped the base interval

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

    private static InvalidDataException DamagedLcp() => new("the LCP array is inconsistent: the file is damaged");
}
