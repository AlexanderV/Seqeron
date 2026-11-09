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
    internal string LongestRepeatedSubstring() => Statistics.LrsDepth == 0 ? string.Empty : Text.Substring(SuffixArray[Statistics.LrsRank], Statistics.LrsDepth);

    /// <summary>The same substring as a slice of the original text: no copy (the index shares the caller's string).</summary>
    internal ReadOnlyMemory<char> LongestRepeatedSubstringMemory() =>
        Statistics.LrsDepth == 0 ? ReadOnlyMemory<char>.Empty : Text.AsMemory(SuffixArray[Statistics.LrsRank], Statistics.LrsDepth);

    /// <summary>
    /// Every distinct longest repeated substring (all length ties), each with all its 0-based start positions ascending (occurrences may overlap), ordered by first occurrence;
    /// empty when no symbol repeats. Identical to the Ukkonen tree's list (<c>SuffixTreeAlgorithms.FindAllLongestRepeatedSubstrings</c>), without the traversal and without a
    /// search per substring: the longest repeats are the lcp-intervals whose value is the maximum LCP, i.e. exactly the maximal runs of that value in the LCP array, and their
    /// positions are the suffix-array slices of the runs. One sequential pass over the LCP array; the positions of all ties together are at most n (the intervals are disjoint);
    /// the substrings are the output (k ties of length L, k &lt;= (n - L + 1) / 2, so at most about n^2 / 8 characters in the worst case, which the return type requires).
    /// </summary>
    internal IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> FindAllLongestRepeatedSubstrings()
    {
        int length = Statistics.LrsDepth;
        if (length == 0) return Array.Empty<(string, IReadOnlyList<int>)>();

        var results = new List<(string Substring, IReadOnlyList<int> Positions)>();
        CompactLcp.Cursor lcp = Lcp.CreateCursor();
        int runStart = -1;      // left boundary of the interval whose run of the maximum value is open
        for (int rank = 1; rank <= Length; rank++)
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
        if (runStart >= 0) results.Add(Tie(runStart, Length, length));

        results.Sort(static (a, b) => a.Positions[0].CompareTo(b.Positions[0]));
        return results;
    }

    private (string Substring, IReadOnlyList<int> Positions) Tie(int lb, int rb, int length)
    {
        int[] positions = SuffixArray.AsSpan(lb, rb - lb + 1).ToArray();
        Array.Sort(positions);
        return (Text.Substring(positions[0], length), positions);
    }

    /// <summary>
    /// The number of distinct non-empty substrings of the text: every suffix contributes its length minus its LCP with the previous suffix in suffix-array order, and the suffix
    /// lengths are 1..n, so the total is <c>n (n + 1) / 2 - sum LCP</c> (in <c>long</c>). One sequential pass over the LCP array, no allocation.
    /// </summary>
    internal long CountDistinctSubstrings()
    {
        long n = Length;
        long total = n * (n + 1) / 2;
        CompactLcp.Cursor lcp = Lcp.CreateCursor();
        for (int rank = 1; rank <= Length; rank++) total -= lcp.Get(rank);
        return total;
    }

    /// <summary>
    /// Element <c>i</c> (<c>1 &lt;= i &lt;= min(maxLength, n)</c>) is the number of distinct substrings of length <c>i</c>, element 0 is 1 (the empty substring; also for the empty
    /// text, which gives <c>[1]</c>). A suffix of length at least <c>i</c> is the first in suffix-array order with its length-<c>i</c> prefix exactly when its LCP with the previous
    /// suffix is below <c>i</c>, and n - i + 1 suffixes have length at least i, so the count is <c>n - i + 1</c> minus the number of ranks whose LCP is at least <c>i</c>: a
    /// histogram of the LCP values (clamped at <c>maxLength</c>) and its suffix sums, built in place in the result array. Two passes, O(n + maxLength) time, no allocation besides
    /// the result (8 bytes per length).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLength"/> is less than 1.</exception>
    internal long[] CountDistinctSubstringsByLength(int maxLength)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLength, 1);
        int n = Length;
        int m = Math.Min(maxLength, n);
        var counts = new long[m + 1];
        CompactLcp.Cursor lcp = Lcp.CreateCursor();
        for (int rank = 1; rank <= n; rank++)
        {
            int value = lcp.Get(rank);
            if (value > 0) counts[Math.Min(value, m)]++;
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
    internal IReadOnlyList<string> GetAllSuffixes()
    {
        var suffixes = new string[Length];
        for (int rank = 1; rank <= Length; rank++) suffixes[rank - 1] = Text.Substring(SuffixArray[rank]);
        return suffixes;
    }

    /// <summary>The suffixes of <see cref="GetAllSuffixes"/> one at a time, lazily: a walk of the suffix array that materialises one suffix per step (nothing is evaluated before the first <c>MoveNext</c>).</summary>
    internal IEnumerable<string> EnumerateSuffixes()
    {
        for (int rank = 1; rank <= Length; rank++) yield return Text.Substring(SuffixArray[rank]);
    }
}
