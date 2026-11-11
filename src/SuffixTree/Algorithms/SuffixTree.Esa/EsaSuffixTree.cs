using SuffixTree.Indexing;

namespace SuffixTree.Esa;

/// <summary>STUB of work package 1b.6 (tests first): every member throws <see cref="NotImplementedException"/>.</summary>
public sealed class EsaSuffixTree : ISuffixTree, ISuffixTreeNavigation, ISuffixArrayAccess, IDisposable
{
    internal EsaSuffixTree(EsaIndex index, EsaStatistics? statistics = null) => throw new NotImplementedException();

    /// <summary>Builds the tree of a text.</summary>
    /// <param name="text">The text.</param>
    /// <returns>The tree.</returns>
    public static EsaSuffixTree Build(string text) => throw new NotImplementedException();

    internal static EsaSuffixTree Build(string text, ChildLayoutKind? layout) => throw new NotImplementedException();

    internal static int CheckedNodeCount(long nodeCount) => throw new NotImplementedException();

    internal static void EnsureTextLength(long textLength) => throw new NotImplementedException();

    internal EsaIndex Index { get; } = null!;

    /// <inheritdoc />
    public ITextSource Text => throw new NotImplementedException();

    /// <inheritdoc />
    public int NodeCount => throw new NotImplementedException();

    /// <inheritdoc />
    public int LeafCount => throw new NotImplementedException();

    /// <inheritdoc />
    public int MaxDepth => throw new NotImplementedException();

    /// <inheritdoc />
    public bool IsEmpty => throw new NotImplementedException();

    /// <inheritdoc />
    public bool Contains(string value) => throw new NotImplementedException();

    /// <inheritdoc />
    public bool Contains(ReadOnlySpan<char> value) => throw new NotImplementedException();

    /// <inheritdoc />
    public IReadOnlyList<int> FindAllOccurrences(string pattern) => throw new NotImplementedException();

    /// <inheritdoc />
    public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern) => throw new NotImplementedException();

    /// <inheritdoc />
    public int CountOccurrences(string pattern) => throw new NotImplementedException();

    /// <inheritdoc />
    public int CountOccurrences(ReadOnlySpan<char> pattern) => throw new NotImplementedException();

    /// <inheritdoc />
    public string LongestRepeatedSubstring() => throw new NotImplementedException();

    /// <inheritdoc />
    public ReadOnlyMemory<char> LongestRepeatedSubstringMemory() => throw new NotImplementedException();

    /// <inheritdoc />
    public IReadOnlyList<string> GetAllSuffixes() => throw new NotImplementedException();

    /// <inheritdoc />
    public IEnumerable<string> EnumerateSuffixes() => throw new NotImplementedException();

    /// <inheritdoc />
    public string LongestCommonSubstring(string other) => throw new NotImplementedException();

    /// <inheritdoc />
    public string LongestCommonSubstring(ReadOnlySpan<char> other) => throw new NotImplementedException();

    /// <inheritdoc />
    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other) => throw new NotImplementedException();

    /// <summary>The same for a character span (a member of the class, not of the interface).</summary>
    /// <param name="other">The other string.</param>
    /// <returns>The longest common substring and its positions.</returns>
    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(ReadOnlySpan<char> other) => Index.LongestCommonSubstringInfo(other);

    /// <inheritdoc />
    public (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other)
        => throw new NotImplementedException();

    /// <inheritdoc />
    public IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)> FindAllDistinctLongestCommonSubstrings(string other)
        => throw new NotImplementedException();

    /// <inheritdoc />
    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(string query, int minLength) => throw new NotImplementedException();

    /// <inheritdoc />
    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindMaximalExactMatches(string query, int minLength) => throw new NotImplementedException();

    /// <inheritdoc />
    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindMaximalUniqueMatches(string query, int minLength, MumUniqueness uniqueness = MumUniqueness.Both)
        => throw new NotImplementedException();

    /// <inheritdoc />
    public long[] CountDistinctSubstringsByLength(int maxLength) => throw new NotImplementedException();

    /// <inheritdoc />
    public long CountDistinctSubstrings() => throw new NotImplementedException();

    /// <inheritdoc />
    public IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> FindAllLongestRepeatedSubstrings() => throw new NotImplementedException();

    /// <inheritdoc />
    public IReadOnlyList<(int FirstPosition, int SecondPosition, int Length)> FindMaximalRepeatedPairs(int minLength, Func<char, bool>? isUniqueSymbol = null)
        => throw new NotImplementedException();

    /// <inheritdoc />
    public string PrintTree() => throw new NotImplementedException();

    /// <inheritdoc />
    public void Traverse(ISuffixTreeVisitor visitor) => throw new NotImplementedException();

    /// <inheritdoc />
    public void Dispose() => throw new NotImplementedException();

    /// <inheritdoc />
    public string Engine => throw new NotImplementedException();

    /// <inheritdoc />
    public long TextLength => throw new NotImplementedException();

    IIndexText ISuffixIndex.Text => throw new NotImplementedException();

    /// <inheritdoc />
    public IndexLocus Root => throw new NotImplementedException();

    /// <inheritdoc />
    public IndexSupport Support(IndexOperation operation) => throw new NotImplementedException();

    /// <inheritdoc />
    public string? SupportNote(IndexOperation operation) => throw new NotImplementedException();

    /// <inheritdoc />
    public long Count(in IndexLocus locus) => throw new NotImplementedException();

    /// <inheritdoc />
    public bool TryExtendRight(in IndexLocus locus, int symbol, out IndexLocus result) => throw new NotImplementedException();

    /// <inheritdoc />
    public bool TryExtendRight(in IndexLocus locus, ReadOnlySpan<char> pattern, out IndexLocus result) => throw new NotImplementedException();

    /// <inheritdoc />
    public bool TryExtendLeft(in IndexLocus locus, int symbol, out IndexLocus result) => throw new NotImplementedException();

    /// <inheritdoc />
    public int Locate(in IndexLocus locus, long skip, Span<long> positions) => throw new NotImplementedException();

    /// <inheritdoc />
    public bool IsExplicitNode(in IndexLocus locus) => throw new NotImplementedException();

    /// <inheritdoc />
    public bool TryParent(in IndexLocus locus, out IndexLocus parent) => throw new NotImplementedException();

    /// <inheritdoc />
    public bool TryChild(in IndexLocus node, int firstSymbol, out IndexLocus child) => throw new NotImplementedException();

    /// <inheritdoc />
    public int GetChildren(in IndexLocus node, Span<IndexLocus> children) => throw new NotImplementedException();

    /// <inheritdoc />
    public bool TrySuffixLink(in IndexLocus node, out IndexLocus link) => throw new NotImplementedException();

    /// <inheritdoc />
    public int EdgeSymbol(in IndexLocus node, long offset) => throw new NotImplementedException();

    /// <inheritdoc />
    public long SuffixAt(long rank) => throw new NotImplementedException();

    /// <inheritdoc />
    public long RankOf(long position) => throw new NotImplementedException();

    /// <inheritdoc />
    public long LcpAt(long rank) => throw new NotImplementedException();

    /// <inheritdoc />
    public long Lcp(long rankA, long rankB) => throw new NotImplementedException();

    /// <inheritdoc />
    public IndexLocus Lca(long rankA, long rankB) => throw new NotImplementedException();
}
