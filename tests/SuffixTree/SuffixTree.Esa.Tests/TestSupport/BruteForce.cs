namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>Definitional pattern search (test oracle L0): all start positions of a pattern by comparing at every position.</summary>
public static class BruteForce
{
    /// <summary>Ascending start positions of <paramref name="pattern"/> in <paramref name="text"/>, overlapping occurrences included; the empty pattern occurs at every position <c>0..n-1</c> (the library convention).</summary>
    public static List<int> Occurrences(string text, string pattern)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(pattern);
        var result = new List<int>();
        if (pattern.Length == 0)
        {
            for (int i = 0; i < text.Length; i++) result.Add(i);
            return result;
        }
        for (int i = 0; i + pattern.Length <= text.Length; i++)
            if (string.CompareOrdinal(text, i, pattern, 0, pattern.Length) == 0) result.Add(i);
        return result;
    }
}
