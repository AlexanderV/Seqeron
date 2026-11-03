namespace SuffixTree.Esa;

/// <summary>
/// The 32-bit instantiation of the reference suffix-array construction and the type/LMS helpers shared by every index width.
/// The algorithm itself is <see cref="SaisReference{TIdx}"/>; this class keeps the short names used by the callers and the tests
/// (<c>SaisReference.BuildSuffixArray</c>, ...) and delegates to <c>SaisReference&lt;int&gt;</c>.
/// </summary>
internal static class SaisReference
{
    /// <summary>Builds the suffix array of a string (ordinal UTF-16 order) through the order-preserving dense remap.</summary>
    internal static int[] BuildSuffixArray(ReadOnlySpan<char> text, SaisTrace? trace = null)
        => SaisReference<int>.BuildSuffixArray(text, trace);

    /// <summary>Builds the suffix array of a symbol sequence over <c>[0, alphabetSize)</c> into a new array of <c>n + 1</c> entries.</summary>
    /// <exception cref="NotSupportedException">The suffix array would not fit the 32-bit path (checked before any allocation).</exception>
    internal static int[] Build<TReader>(TReader text, int alphabetSize, SaisTrace? trace = null)
        where TReader : ISymbolReader, allows ref struct
        => SaisReference<int>.Build(text, alphabetSize, trace);

    /// <summary>Builds the suffix array into <paramref name="sa"/>, which must have exactly <c>n + 1</c> entries.</summary>
    internal static void Build<TReader>(TReader text, int alphabetSize, Span<int> sa, SaisTrace? trace = null)
        where TReader : ISymbolReader, allows ref struct
        => SaisReference<int>.Build(text, alphabetSize, sa, trace);

    /// <summary>Type of every position, <c>true</c> = S-type. The virtual sentinel at index n is S-type, so position n - 1 is always L-type.</summary>
    internal static bool[] ClassifyTypes<TReader>(TReader t)
        where TReader : ISymbolReader, allows ref struct
    {
        int n = t.Length;
        var isS = new bool[n + 1];
        isS[n] = true;
        for (int i = n - 2; i >= 0; i--)
        {
            int a = t[i], b = t[i + 1];
            isS[i] = a < b || (a == b && isS[i + 1]);
        }
        return isS;
    }

    /// <summary>True when position i is an LMS position: S-type with an L-type predecessor (position 0 never is).</summary>
    internal static bool IsLms(bool[] isS, int i) => i > 0 && isS[i] && !isS[i - 1];

    /// <summary>Number of LMS positions in <c>[1, n - 1]</c>; the sentinel position n is the (implicit) smallest LMS position and is not counted.</summary>
    internal static int CountLms(bool[] isS)
    {
        int m = 0;
        for (int i = 1; i < isS.Length - 1; i++)
            if (IsLms(isS, i)) m++;
        return m;
    }

    /// <summary>
    /// Equality of the LMS substrings starting at the distinct LMS positions a and b (same length, same symbols, same types).
    /// The LMS substring of a is text[a..next LMS position] inclusive; the one that reaches the sentinel is unique.
    /// Internal so that the tests can check it against the definition on its own.
    /// </summary>
    internal static bool LmsSubstringsEqual<TReader>(TReader t, bool[] isS, int a, int b)
        where TReader : ISymbolReader, allows ref struct
    {
        int n = t.Length;
        for (int d = 0; a + d <= n && b + d <= n; d++)
        {
            int x = a + d, y = b + d;
            if (x == n || y == n) return false; // a != b, so at most one of them is the sentinel, which no other substring contains
            if (t[x] != t[y] || isS[x] != isS[y]) return false;
            if (d > 0 && IsLms(isS, x)) return true; // equal so far, so y is LMS as well: both substrings end here
        }
        return false; // not reached: the loop ends at the sentinel of the later substring
    }
}
