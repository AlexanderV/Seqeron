namespace SuffixTree.Esa;

/// <summary>
/// Native analysis operations of the in-memory ESA (work package 1b.3): the longest repeated substring, all of its ties, the distinct-substring counts and the suffix
/// enumeration, computed from the suffix array, the LCP array and the header statistics alone. Signatures and results are those of <c>ISuffixTreeAnalysis</c> and
/// <c>ISuffixTree</c> (the facade of work package 1b.6 forwards to them); the maximal repeated pairs and the k-common-substring operations are NOT here: they run on
/// <see cref="Traverse"/> through the generic algorithms of <c>SuffixTree.Core</c>.
/// </summary>
internal sealed partial class EsaIndex
{
    /// <summary>
    /// A longest substring that occurs at least twice (occurrences may overlap), or the empty string when no symbol repeats. It is the path label of a deepest internal node,
    /// which the construction scan recorded in the header (<see cref="EsaStatistics.LrsDepth"/>, <see cref="EsaStatistics.LrsRank"/>): O(L) to copy, no search. On a length tie
    /// the contract leaves the choice open (<c>ISuffixTreeAnalysis.cs:13</c>); the ESA returns the lexicographically smallest of the tied substrings (the deepest interval that
    /// closes first, <see cref="StatisticsSink"/>). <see cref="FindAllLongestRepeatedSubstrings"/> lists every tie.
    /// </summary>
    /// <exception cref="NotSupportedException">The longest repeat has more than <see cref="int.MaxValue"/> symbols (it is no string): <see cref="LrsLength"/> and <see cref="LrsPosition"/> are the 64-bit members.</exception>
    internal string LongestRepeatedSubstring() =>
        Statistics.LrsDepth == 0 ? string.Empty : TextView.Substring(Sa[Statistics.LrsRank], IntBridge.Length(Statistics.LrsDepth, "EsaIndex.LongestRepeatedSubstring (the string of the longest repeat; the 64-bit members are LrsLength and LrsPosition)"));

    /// <summary>The same substring as a slice of the original text: no copy (the index shares the caller's string).</summary>
    /// <exception cref="NotSupportedException">The longest repeat has more than <see cref="int.MaxValue"/> symbols.</exception>
    internal ReadOnlyMemory<char> LongestRepeatedSubstringMemory() =>
        Statistics.LrsDepth == 0 ? ReadOnlyMemory<char>.Empty : TextView.AsMemory(Sa[Statistics.LrsRank], IntBridge.Length(Statistics.LrsDepth, "EsaIndex.LongestRepeatedSubstringMemory (the memory of the longest repeat; the 64-bit members are LrsLength and LrsPosition)"));

    /// <summary>The length of the longest repeated substring (0 when no symbol repeats), at 64 bits: a text above 2^31 symbols has repeats above 2^31 (<c>a^n</c>: n - 1).</summary>
    internal long LrsLength => Statistics.LrsDepth;

    /// <summary>The start position of the longest repeated substring (<see cref="LongestRepeatedSubstring"/>) at 64 bits, -1 when no symbol repeats.</summary>
    internal long LrsPosition => Statistics.LrsDepth == 0 ? -1 : Sa[Statistics.LrsRank];

    /// <summary>
    /// Every distinct longest repeated substring (all length ties), each with all its 0-based start positions ascending (occurrences may overlap), ordered by first occurrence;
    /// empty when no symbol repeats. Identical to the Ukkonen tree's list (<c>SuffixTreeAlgorithms.FindAllLongestRepeatedSubstrings</c>), without the traversal and without a
    /// search per substring: the longest repeats are the lcp-intervals whose value is the maximum LCP, i.e. exactly the maximal runs of that value in the LCP array, and their
    /// positions are the suffix-array slices of the runs. One sequential pass over the LCP array; the positions of all ties together are at most n (the intervals are disjoint);
    /// the substrings are the output (k ties of length L, k &lt;= (n - L + 1) / 2, so at most about n^2 / 8 characters in the worst case, which the return type requires).
    /// </summary>
    /// <exception cref="NotSupportedException">A repeat, or the list of its positions, does not fit the 32-bit shape of the result.</exception>
    internal IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> FindAllLongestRepeatedSubstrings()
    {
        long length = Statistics.LrsDepth;
        if (length == 0) return Array.Empty<(string, IReadOnlyList<int>)>();

        var results = new List<(string Substring, IReadOnlyList<int> Positions)>();
        CompactLcp.LongCursor lcp = Lcp.CreateLongCursor();
        long runStart = -1;      // left boundary of the interval whose run of the maximum value is open
        for (long rank = 1; rank <= _length; rank++)
        {
            if (lcp.Get(rank) == length)
            {
                if (runStart < 0) runStart = rank - 1;
            }
            else if (runStart >= 0)
            {
                results.Add(Tie(runStart, rank - 1, length));
                runStart = -1;
            }
        }
        if (runStart >= 0) results.Add(Tie(runStart, _length, length));

        results.Sort(static (a, b) => a.Positions[0].CompareTo(b.Positions[0]));
        return results;
    }

    private (string Substring, IReadOnlyList<int> Positions) Tie(long lb, long rb, long length)
    {
        int[] positions = Sa.ToArray(lb, IntBridge.Length(rb - lb + 1, "EsaIndex.FindAllLongestRepeatedSubstrings (the positions of one longest repeat)"));
        Array.Sort(positions);
        return (TextView.Substring(positions[0], IntBridge.Length(length, "EsaIndex.FindAllLongestRepeatedSubstrings (the string of a longest repeat)")), positions);
    }

    /// <summary>
    /// The number of distinct non-empty substrings of the text: every suffix contributes its length minus its LCP with the previous suffix in suffix-array order, and the suffix
    /// lengths are 1..n, so the total is <c>n (n + 1) / 2 - sum LCP</c> (in 128 bits: the two terms are of the order of n^2). One sequential pass over the LCP array, no allocation.
    /// </summary>
    /// <exception cref="NotSupportedException">The number does not fit 64 bits (a text of more than about 2^32 symbols with few repeats).</exception>
    internal long CountDistinctSubstrings() => CountDistinctSubstrings(Lcp, _length);

    /// <summary>The distinct substrings of a text of <paramref name="n"/> symbols whose LCP array is <paramref name="lcp"/> (any size and kind): <see cref="CountDistinctSubstrings()"/> without an index around it.</summary>
    /// <exception cref="NotSupportedException">The number does not fit 64 bits.</exception>
    internal static long CountDistinctSubstrings(CompactLcp lcp, long n)
    {
        Int128 length = n;
        Int128 total = length * (length + 1) / 2;
        CompactLcp.LongCursor cursor = lcp.CreateLongCursor();
        for (long rank = 1; rank <= n; rank++) total -= cursor.Get(rank);
        if (total > long.MaxValue) throw new NotSupportedException($"the number of distinct substrings is {total}, more than 64 bits hold");
        return (long)total;
    }

    /// <summary>
    /// Element <c>i</c> (<c>1 &lt;= i &lt;= min(maxLength, n)</c>) is the number of distinct substrings of length <c>i</c>, element 0 is 1 (the empty substring; also for the empty
    /// text, which gives <c>[1]</c>). A suffix of length at least <c>i</c> is the first in suffix-array order with its length-<c>i</c> prefix exactly when its LCP with the previous
    /// suffix is below <c>i</c>, and n - i + 1 suffixes have length at least i, so the count is <c>n - i + 1</c> minus the number of ranks whose LCP is at least <c>i</c>: a
    /// histogram of the LCP values (clamped at <c>maxLength</c>) and its suffix sums, built in place in the result array. Two passes, O(n + maxLength) time, no allocation besides
    /// the result (8 bytes per length).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLength"/> is less than 1.</exception>
    internal long[] CountDistinctSubstringsByLength(int maxLength) => CountDistinctSubstringsByLength(Lcp, _length, maxLength);

    /// <summary>The same for an LCP array of any size and kind, without an index around it.</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLength"/> is less than 1.</exception>
    internal static long[] CountDistinctSubstringsByLength(CompactLcp lcp, long n, int maxLength)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLength, 1);
        int m = n < maxLength ? (int)n : maxLength; // (int): below maxLength, which is an int
        var counts = new long[m + 1];
        CompactLcp.LongCursor cursor = lcp.CreateLongCursor();
        for (long rank = 1; rank <= n; rank++)
        {
            long value = cursor.Get(rank);
            if (value > 0) counts[value < m ? (int)value : m]++; // (int): below m, which is an int
        }

        long atLeast = 0;       // number of ranks whose LCP is at least i
        for (int i = m; i >= 1; i--)
        {
            atLeast += counts[i];
            counts[i] = n - i + 1 - atLeast;
        }
        counts[0] = 1;
        return counts;
    }

    /// <summary>
    /// All suffixes of the text in lexicographic (ordinal) order, the order of the suffix array: a prefix sorts before its extensions, as in the Ukkonen tree (terminator first).
    /// The empty suffix is not listed. The result holds n strings of average length n / 2, i.e. Theta(n^2) characters by definition of the operation; use
    /// <see cref="EnumerateSuffixes"/> to avoid materialising them all.
    /// </summary>
    /// <exception cref="NotSupportedException">The text has more than <see cref="int.MaxValue"/> symbols (n strings do not fit an array; <see cref="EnumerateSuffixes"/> is lazy), or a suffix is longer than a string.</exception>
    internal IReadOnlyList<string> GetAllSuffixes()
    {
        var suffixes = new string[IntBridge.Length(_length, "EsaIndex.GetAllSuffixes (the list holds n strings; EnumerateSuffixes is lazy)")];
        SaView sa = Sa;
        for (int rank = 1; rank <= suffixes.Length; rank++) suffixes[rank - 1] = TextView.Substring(sa[(long)rank]);
        return suffixes;
    }

    /// <summary>The suffixes of <see cref="GetAllSuffixes"/> one at a time, lazily: a walk of the suffix array that materialises one suffix per step (nothing is evaluated before the first <c>MoveNext</c>).</summary>
    /// <exception cref="NotSupportedException">A suffix is longer than a string (a text above 2^31 symbols).</exception>
    internal IEnumerable<string> EnumerateSuffixes()
    {
        SaView sa = Sa;
        for (long rank = 1; rank <= _length; rank++) yield return TextView.Substring(sa[rank]);
    }
}
