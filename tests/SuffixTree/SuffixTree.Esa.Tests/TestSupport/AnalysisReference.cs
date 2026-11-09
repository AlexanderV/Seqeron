namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>
/// Definition-level references (test oracle L0) of the analysis operations of work package 1b.3. They are written from the definitions only (sets of substrings,
/// comparison of every pair), share no code with the library, and are meant for texts of at most a few hundred symbols.
/// </summary>
public static class AnalysisReference
{
    /// <summary>
    /// Element <c>i</c> (<c>1 &lt;= i &lt;= min(maxLength, n)</c>) is the number of distinct substrings of length <c>i</c>; element 0 is 1 (the empty substring), also for the empty text.
    /// </summary>
    public static long[] DistinctByLength(string text, int maxLength)
    {
        int m = Math.Min(maxLength, text.Length);
        var counts = new long[m + 1];
        counts[0] = 1;
        for (int length = 1; length <= m; length++)
        {
            var seen = new HashSet<string>(StringComparer.Ordinal);
            for (int start = 0; start + length <= text.Length; start++) seen.Add(text.Substring(start, length));
            counts[length] = seen.Count;
        }
        return counts;
    }

    /// <summary>The number of distinct non-empty substrings.</summary>
    public static long DistinctTotal(string text)
    {
        var seen = new HashSet<string>(StringComparer.Ordinal);
        for (int start = 0; start < text.Length; start++)
            for (int length = 1; start + length <= text.Length; length++) seen.Add(text.Substring(start, length));
        return seen.Count;
    }

    /// <summary>
    /// Every distinct substring of the greatest length that occurs at least twice (overlapping occurrences count), with all its start positions ascending, ordered by the
    /// first occurrence. Empty when no symbol repeats.
    /// </summary>
    public static List<(string Substring, int[] Positions)> AllLongestRepeated(string text)
    {
        var result = new List<(string Substring, int[] Positions)>();
        for (int length = text.Length - 1; length >= 1 && result.Count == 0; length--)
        {
            var positions = new Dictionary<string, List<int>>(StringComparer.Ordinal);
            for (int start = 0; start + length <= text.Length; start++)
            {
                string s = text.Substring(start, length);
                if (!positions.TryGetValue(s, out List<int>? list)) positions[s] = list = [];
                list.Add(start);
            }
            foreach (var (s, list) in positions)
                if (list.Count >= 2) result.Add((s, list.ToArray()));
        }
        result.Sort((a, b) => a.Positions[0].CompareTo(b.Positions[0]));
        return result;
    }

    /// <summary>All suffixes of the text, non-empty, in ordinal (UTF-16 code unit) order; a prefix sorts before its extensions.</summary>
    public static string[] SortedSuffixes(string text)
    {
        var suffixes = new string[text.Length];
        for (int i = 0; i < suffixes.Length; i++) suffixes[i] = text.Substring(i);
        Array.Sort(suffixes, StringComparer.Ordinal);
        return suffixes;
    }

    /// <summary>
    /// Maximal repeated pairs by the definition: every <c>i &lt; j</c> whose match, extended to the right as far as the texts agree, has length at least <paramref name="minLength"/> and
    /// that cannot be extended to the left. A symbol for which <paramref name="isUniqueSymbol"/> returns true matches nothing, not even itself. Ordered by <c>i</c>, then <c>j</c>.
    /// </summary>
    public static List<(int FirstPosition, int SecondPosition, int Length)> MaximalRepeatedPairs(string text, int minLength, Func<char, bool>? isUniqueSymbol)
    {
        bool Same(int a, int b) => text[a] == text[b] && (isUniqueSymbol is null || !isUniqueSymbol(text[a]));
        var pairs = new List<(int, int, int)>();
        for (int i = 0; i < text.Length; i++)
            for (int j = i + 1; j < text.Length; j++)
            {
                if (i > 0 && Same(i - 1, j - 1)) continue;          // extendable to the left
                int length = 0;
                while (j + length < text.Length && Same(i + length, j + length)) length++;
                if (length >= minLength) pairs.Add((i, j, length));
            }
        return pairs;
    }

    /// <summary>
    /// Gusfield's l(q) for q = 1..k as element q (element 0 equals element 1): the length of the longest substring that occurs in at least q of the texts, and for
    /// each q the distinct substrings of that length with that support, ordinally sorted (empty when the length is 0).
    /// </summary>
    public static (int[] Lengths, string[][] BySupport) CommonSubstrings(IReadOnlyList<string> texts)
    {
        int k = texts.Count;
        int longest = texts.Max(t => t.Length);
        var support = new Dictionary<string, int>(StringComparer.Ordinal);
        var lengths = new int[k + 1];
        var bySupport = new string[k + 1][];
        for (int q = 0; q <= k; q++) bySupport[q] = [];
        for (int length = longest; length >= 1; length--)
        {
            support.Clear();
            foreach (string t in texts)
                foreach (string s in Enumerable.Range(0, Math.Max(0, t.Length - length + 1)).Select(i => t.Substring(i, length)).Distinct(StringComparer.Ordinal))
                    support[s] = support.GetValueOrDefault(s) + 1;
            for (int q = 1; q <= k; q++)
            {
                if (lengths[q] != 0) continue;
                string[] hits = support.Where(kv => kv.Value >= q).Select(kv => kv.Key).Order(StringComparer.Ordinal).ToArray();
                if (hits.Length == 0) continue;
                lengths[q] = length;
                bySupport[q] = hits;
            }
        }
        lengths[0] = lengths[1];
        bySupport[0] = bySupport[1];
        return (lengths, bySupport);
    }
}
