using SuffixTree.Esa;

namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>
/// Test-only <see cref="ISuffixTree"/> view of an <see cref="EsaIndex"/> with just what the oracles need: the text, the statistics and
/// <see cref="Traverse"/>, so that <c>SuffixTreeSerializer.CalculateLogicalHash</c> (which takes an <see cref="ISuffixTree"/>) can be applied to the
/// ESA. The facade of the product (work package 1b.6) replaces it; every member outside that set throws <see cref="NotSupportedException"/>.
/// </summary>
internal sealed class EsaTreeShim(EsaIndex index) : ISuffixTree
{
    public ITextSource Text { get; } = new StringTextSource(index.Text);

    public int NodeCount => checked((int)index.Statistics.NodeCount);

    public int LeafCount => index.Statistics.LeafCount;

    public int MaxDepth => index.Statistics.MaxDepth;

    public bool IsEmpty => index.Length == 0;

    public void Traverse(ISuffixTreeVisitor visitor) => index.Traverse(visitor);

    public string PrintTree() => throw Unsupported();

    public bool Contains(string value) => throw Unsupported();

    public bool Contains(ReadOnlySpan<char> value) => throw Unsupported();

    public IReadOnlyList<int> FindAllOccurrences(string pattern) => throw Unsupported();

    public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern) => throw Unsupported();

    public int CountOccurrences(string pattern) => throw Unsupported();

    public int CountOccurrences(ReadOnlySpan<char> pattern) => throw Unsupported();

    public string LongestRepeatedSubstring() => throw Unsupported();

    public IReadOnlyList<string> GetAllSuffixes() => throw Unsupported();

    public IEnumerable<string> EnumerateSuffixes() => throw Unsupported();

    public string LongestCommonSubstring(string other) => throw Unsupported();

    public string LongestCommonSubstring(ReadOnlySpan<char> other) => throw Unsupported();

    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other) => throw Unsupported();

    public (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other)
        => throw Unsupported();

    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(string query, int minLength) => throw Unsupported();

    private static NotSupportedException Unsupported() => new("the test shim of the ESA implements only Text, the statistics and Traverse");
}
