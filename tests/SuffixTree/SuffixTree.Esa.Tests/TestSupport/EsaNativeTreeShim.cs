using SuffixTree.Esa;

namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>
/// Test-only <see cref="ISuffixTree"/> view of an <see cref="EsaIndex"/> as the facade of work package 1b.6 will be: the text, the statistics, <see cref="Traverse"/>, the search members and the
/// NATIVE analysis operations of work package 1b.3 are forwarded to the index. The members that the interface implements by default through <c>Traverse</c>
/// (<c>FindMaximalRepeatedPairs</c>) are deliberately not declared, so a call through <see cref="ISuffixTree"/> takes the generic path; <see cref="EsaTreeShim"/> is the view
/// that takes the generic path for ALL of them (the second oracle of the native overrides). The matching-statistics members (1b.4) throw <see cref="NotSupportedException"/>.
/// </summary>
internal sealed class EsaNativeTreeShim(EsaIndex index) : ISuffixTree
{
    public ITextSource Text { get; } = new StringTextSource(index.Text);

    public int NodeCount => checked((int)index.Statistics.NodeCount);

    public int LeafCount => index.Statistics.LeafCount;

    public int MaxDepth => index.Statistics.MaxDepth;

    public bool IsEmpty => index.Length == 0;

    public void Traverse(ISuffixTreeVisitor visitor) => index.Traverse(visitor);

    public string PrintTree() => throw Unsupported();

    public bool Contains(string value) => index.Contains(value);

    public bool Contains(ReadOnlySpan<char> value) => index.Contains(value);

    public IReadOnlyList<int> FindAllOccurrences(string pattern) => index.FindAllOccurrences(pattern);

    public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern) => index.FindAllOccurrences(pattern);

    public int CountOccurrences(string pattern) => index.CountOccurrences(pattern);

    public int CountOccurrences(ReadOnlySpan<char> pattern) => index.CountOccurrences(pattern);

    public string LongestRepeatedSubstring() => index.LongestRepeatedSubstring();

    public ReadOnlyMemory<char> LongestRepeatedSubstringMemory() => index.LongestRepeatedSubstringMemory();

    public IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> FindAllLongestRepeatedSubstrings() => index.FindAllLongestRepeatedSubstrings();

    public long CountDistinctSubstrings() => index.CountDistinctSubstrings();

    public long[] CountDistinctSubstringsByLength(int maxLength) => index.CountDistinctSubstringsByLength(maxLength);

    public IReadOnlyList<string> GetAllSuffixes() => index.GetAllSuffixes();

    public IEnumerable<string> EnumerateSuffixes() => index.EnumerateSuffixes();

    public string LongestCommonSubstring(string other) => throw Unsupported();

    public string LongestCommonSubstring(ReadOnlySpan<char> other) => throw Unsupported();

    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other) => throw Unsupported();

    public (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other)
        => throw Unsupported();

    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(string query, int minLength) => throw Unsupported();

    private static NotSupportedException Unsupported() => new("the test shim implements the search members, Traverse and the native analysis operations of 1b.3 only");
}
