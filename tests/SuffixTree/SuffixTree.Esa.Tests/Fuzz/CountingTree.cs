using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// A decorator over <see cref="EsaSuffixTree"/> that forwards every member of the audited interfaces and counts it under its key (<see cref="FuzzMemberKeys"/>) in a <see cref="MemberCounters"/>: the means by which the
/// harness reports per-member call counts and by which a test proves that every member is called on the in-memory and on the file-backed tree. The members with a default implementation in <see cref="ISuffixTree"/> are
/// forwarded explicitly (an inherited default would silently skip the engine's own override), what a member returns (<see cref="ITextSource"/>, <see cref="IIndexText"/>, the visitor of <c>Traverse</c>) is wrapped so that its
/// members are counted as well.
/// </summary>
internal sealed class CountingTree : ISuffixTree, ISuffixTreeNavigation, ISuffixArrayAccess, IDisposable
{
    private readonly EsaSuffixTree _inner;
    private readonly MemberCounters _counters;
    private readonly int _engine;
    private readonly CountingText _text;
    private readonly CountingIndexText _indexText;

    internal CountingTree(EsaSuffixTree inner, MemberCounters counters, int engine)
    {
        _inner = inner;
        _counters = counters;
        _engine = engine;
        _text = new CountingText(this);
        _indexText = new CountingIndexText(this);
    }

    /// <summary>The tree that is counted.</summary>
    internal EsaSuffixTree Inner => _inner;

    /// <summary>Fault injection of the harness tests: a function from the pattern and the true count to the count that this tree reports.</summary>
    internal Func<string, int, int>? CountFault { get; set; }

    internal void Hit(string key) => _counters.Hit(_engine, key);

    // ---------------------------------------------------------------- ISuffixTreeSearch

    public ITextSource Text
    {
        get
        {
            Hit("ISuffixTreeSearch.Text");
            _ = _inner.Text;
            return _text;
        }
    }

    public int NodeCount
    {
        get
        {
            Hit("ISuffixTreeSearch.NodeCount");
            return _inner.NodeCount;
        }
    }

    public int LeafCount
    {
        get
        {
            Hit("ISuffixTreeSearch.LeafCount");
            return _inner.LeafCount;
        }
    }

    public int MaxDepth
    {
        get
        {
            Hit("ISuffixTreeSearch.MaxDepth");
            return _inner.MaxDepth;
        }
    }

    public bool IsEmpty
    {
        get
        {
            Hit("ISuffixTreeSearch.IsEmpty");
            return _inner.IsEmpty;
        }
    }

    public bool Contains(string value)
    {
        Hit("ISuffixTreeSearch.Contains(String)");
        return _inner.Contains(value);
    }

    public bool Contains(ReadOnlySpan<char> value)
    {
        Hit("ISuffixTreeSearch.Contains(ReadOnlySpan<Char>)");
        return _inner.Contains(value);
    }

    public IReadOnlyList<int> FindAllOccurrences(string pattern)
    {
        Hit("ISuffixTreeSearch.FindAllOccurrences(String)");
        return _inner.FindAllOccurrences(pattern);
    }

    public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern)
    {
        Hit("ISuffixTreeSearch.FindAllOccurrences(ReadOnlySpan<Char>)");
        return _inner.FindAllOccurrences(pattern);
    }

    public int CountOccurrences(string pattern)
    {
        Hit("ISuffixTreeSearch.CountOccurrences(String)");
        int count = _inner.CountOccurrences(pattern);
        return CountFault is null ? count : CountFault(pattern, count);
    }

    public int CountOccurrences(ReadOnlySpan<char> pattern)
    {
        Hit("ISuffixTreeSearch.CountOccurrences(ReadOnlySpan<Char>)");
        return _inner.CountOccurrences(pattern);
    }

    // ---------------------------------------------------------------- ISuffixTreeAnalysis

    public string LongestRepeatedSubstring()
    {
        Hit("ISuffixTreeAnalysis.LongestRepeatedSubstring()");
        return _inner.LongestRepeatedSubstring();
    }

    public ReadOnlyMemory<char> LongestRepeatedSubstringMemory()
    {
        Hit("ISuffixTreeAnalysis.LongestRepeatedSubstringMemory()");
        return _inner.LongestRepeatedSubstringMemory();
    }

    public IReadOnlyList<string> GetAllSuffixes()
    {
        Hit("ISuffixTreeAnalysis.GetAllSuffixes()");
        return _inner.GetAllSuffixes();
    }

    public IEnumerable<string> EnumerateSuffixes()
    {
        Hit("ISuffixTreeAnalysis.EnumerateSuffixes()");
        return _inner.EnumerateSuffixes();
    }

    public string LongestCommonSubstring(string other)
    {
        Hit("ISuffixTreeAnalysis.LongestCommonSubstring(String)");
        return _inner.LongestCommonSubstring(other);
    }

    public string LongestCommonSubstring(ReadOnlySpan<char> other)
    {
        Hit("ISuffixTreeAnalysis.LongestCommonSubstring(ReadOnlySpan<Char>)");
        return _inner.LongestCommonSubstring(other);
    }

    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other)
    {
        Hit("ISuffixTreeAnalysis.LongestCommonSubstringInfo(String)");
        return _inner.LongestCommonSubstringInfo(other);
    }

    public (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other)
    {
        Hit("ISuffixTreeAnalysis.FindAllLongestCommonSubstrings(String)");
        return _inner.FindAllLongestCommonSubstrings(other);
    }

    public IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)> FindAllDistinctLongestCommonSubstrings(string other)
    {
        Hit("ISuffixTreeAnalysis.FindAllDistinctLongestCommonSubstrings(String)");
        return _inner.FindAllDistinctLongestCommonSubstrings(other);
    }

    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(string query, int minLength)
    {
        Hit("ISuffixTreeAnalysis.FindExactMatchAnchors(String,Int32)");
        return _inner.FindExactMatchAnchors(query, minLength);
    }

    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindMaximalExactMatches(string query, int minLength)
    {
        Hit("ISuffixTreeAnalysis.FindMaximalExactMatches(String,Int32)");
        return _inner.FindMaximalExactMatches(query, minLength);
    }

    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindMaximalUniqueMatches(string query, int minLength, MumUniqueness uniqueness = MumUniqueness.Both)
    {
        Hit("ISuffixTreeAnalysis.FindMaximalUniqueMatches(String,Int32,MumUniqueness)");
        return _inner.FindMaximalUniqueMatches(query, minLength, uniqueness);
    }

    // ---------------------------------------------------------------- ISuffixTree

    public long[] CountDistinctSubstringsByLength(int maxLength)
    {
        Hit("ISuffixTree.CountDistinctSubstringsByLength(Int32)");
        return _inner.CountDistinctSubstringsByLength(maxLength);
    }

    public long CountDistinctSubstrings()
    {
        Hit("ISuffixTree.CountDistinctSubstrings()");
        return _inner.CountDistinctSubstrings();
    }

    public IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> FindAllLongestRepeatedSubstrings()
    {
        Hit("ISuffixTree.FindAllLongestRepeatedSubstrings()");
        return _inner.FindAllLongestRepeatedSubstrings();
    }

    public IReadOnlyList<(int FirstPosition, int SecondPosition, int Length)> FindMaximalRepeatedPairs(int minLength, Func<char, bool>? isUniqueSymbol = null)
    {
        Hit("ISuffixTree.FindMaximalRepeatedPairs(Int32,Func<Char,Boolean>)");
        return _inner.FindMaximalRepeatedPairs(minLength, isUniqueSymbol);
    }

    // ---------------------------------------------------------------- ISuffixTreeDiagnostics

    public string PrintTree()
    {
        Hit("ISuffixTreeDiagnostics.PrintTree()");
        return _inner.PrintTree();
    }

    public void Traverse(ISuffixTreeVisitor visitor)
    {
        Hit("ISuffixTreeDiagnostics.Traverse(ISuffixTreeVisitor)");
        _inner.Traverse(visitor is null ? null! : new CountingVisitor(this, visitor)); // a null visitor reaches the tree itself, which must refuse it
    }

    // ---------------------------------------------------------------- ISuffixIndex

    public string Engine
    {
        get
        {
            Hit("ISuffixIndex.Engine");
            return _inner.Engine;
        }
    }

    public long TextLength
    {
        get
        {
            Hit("ISuffixIndex.TextLength");
            return _inner.TextLength;
        }
    }

    IIndexText ISuffixIndex.Text
    {
        get
        {
            Hit("ISuffixIndex.Text");
            _ = ((ISuffixIndex)_inner).Text;
            return _indexText;
        }
    }

    public IndexLocus Root
    {
        get
        {
            Hit("ISuffixIndex.Root");
            return _inner.Root;
        }
    }

    public IndexSupport Support(IndexOperation operation)
    {
        Hit("ISuffixIndex.Support(IndexOperation)");
        return _inner.Support(operation);
    }

    public string? SupportNote(IndexOperation operation)
    {
        Hit("ISuffixIndex.SupportNote(IndexOperation)");
        return _inner.SupportNote(operation);
    }

    public long Count(in IndexLocus locus)
    {
        Hit("ISuffixIndex.Count(IndexLocus)");
        return _inner.Count(in locus);
    }

    public bool TryExtendRight(in IndexLocus locus, int symbol, out IndexLocus result)
    {
        Hit("ISuffixIndex.TryExtendRight(IndexLocus,Int32,IndexLocus)");
        return _inner.TryExtendRight(in locus, symbol, out result);
    }

    public bool TryExtendRight(in IndexLocus locus, ReadOnlySpan<char> pattern, out IndexLocus result)
    {
        Hit("ISuffixIndex.TryExtendRight(IndexLocus,ReadOnlySpan<Char>,IndexLocus)");
        return _inner.TryExtendRight(in locus, pattern, out result);
    }

    public bool TryExtendLeft(in IndexLocus locus, int symbol, out IndexLocus result)
    {
        Hit("ISuffixIndex.TryExtendLeft(IndexLocus,Int32,IndexLocus)");
        return _inner.TryExtendLeft(in locus, symbol, out result);
    }

    public int Locate(in IndexLocus locus, long skip, Span<long> positions)
    {
        Hit("ISuffixIndex.Locate(IndexLocus,Int64,Span<Int64>)");
        return _inner.Locate(in locus, skip, positions);
    }

    // ---------------------------------------------------------------- ISuffixTreeNavigation

    public bool IsExplicitNode(in IndexLocus locus)
    {
        Hit("ISuffixTreeNavigation.IsExplicitNode(IndexLocus)");
        return _inner.IsExplicitNode(in locus);
    }

    public bool TryParent(in IndexLocus locus, out IndexLocus parent)
    {
        Hit("ISuffixTreeNavigation.TryParent(IndexLocus,IndexLocus)");
        return _inner.TryParent(in locus, out parent);
    }

    public bool TryChild(in IndexLocus node, int firstSymbol, out IndexLocus child)
    {
        Hit("ISuffixTreeNavigation.TryChild(IndexLocus,Int32,IndexLocus)");
        return _inner.TryChild(in node, firstSymbol, out child);
    }

    public int GetChildren(in IndexLocus node, Span<IndexLocus> children)
    {
        Hit("ISuffixTreeNavigation.GetChildren(IndexLocus,Span<IndexLocus>)");
        return _inner.GetChildren(in node, children);
    }

    public bool TrySuffixLink(in IndexLocus node, out IndexLocus link)
    {
        Hit("ISuffixTreeNavigation.TrySuffixLink(IndexLocus,IndexLocus)");
        return _inner.TrySuffixLink(in node, out link);
    }

    public int EdgeSymbol(in IndexLocus node, long offset)
    {
        Hit("ISuffixTreeNavigation.EdgeSymbol(IndexLocus,Int64)");
        return _inner.EdgeSymbol(in node, offset);
    }

    // ---------------------------------------------------------------- ISuffixArrayAccess

    public long SuffixAt(long rank)
    {
        Hit("ISuffixArrayAccess.SuffixAt(Int64)");
        return _inner.SuffixAt(rank);
    }

    public long RankOf(long position)
    {
        Hit("ISuffixArrayAccess.RankOf(Int64)");
        return _inner.RankOf(position);
    }

    public long LcpAt(long rank)
    {
        Hit("ISuffixArrayAccess.LcpAt(Int64)");
        return _inner.LcpAt(rank);
    }

    public long Lcp(long rankA, long rankB)
    {
        Hit("ISuffixArrayAccess.Lcp(Int64,Int64)");
        return _inner.Lcp(rankA, rankB);
    }

    public IndexLocus Lca(long rankA, long rankB)
    {
        Hit("ISuffixArrayAccess.Lca(Int64,Int64)");
        return _inner.Lca(rankA, rankB);
    }

    // ---------------------------------------------------------------- IDisposable

    public void Dispose()
    {
        Hit("IDisposable.Dispose()");
        _inner.Dispose();
    }

    private sealed class CountingText(CountingTree owner) : ITextSource
    {
        private ITextSource Source => owner._inner.Text;

        public int Length
        {
            get
            {
                owner.Hit("ITextSource.Length");
                return Source.Length;
            }
        }

        public char this[int index]
        {
            get
            {
                owner.Hit("ITextSource.Item[Int32]");
                return Source[index];
            }
        }

        public string Substring(int start, int length)
        {
            owner.Hit("ITextSource.Substring(Int32,Int32)");
            return Source.Substring(start, length);
        }

        public ReadOnlySpan<char> Slice(int start, int length)
        {
            owner.Hit("ITextSource.Slice(Int32,Int32)");
            return Source.Slice(start, length);
        }

        public IEnumerator<char> GetEnumerator() => Source.GetEnumerator();

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() => GetEnumerator();

        public override string ToString() => Source.ToString() ?? string.Empty;
    }

    private sealed class CountingIndexText(CountingTree owner) : IIndexText
    {
        private IIndexText Source => ((ISuffixIndex)owner._inner).Text;

        public long Length
        {
            get
            {
                owner.Hit("IIndexText.Length");
                return Source.Length;
            }
        }

        public int SymbolAt(long position)
        {
            owner.Hit("IIndexText.SymbolAt(Int64)");
            return Source.SymbolAt(position);
        }

        public int CopyTo(long start, Span<char> destination)
        {
            owner.Hit("IIndexText.CopyTo(Int64,Span<Char>)");
            return Source.CopyTo(start, destination);
        }
    }

    private sealed class CountingVisitor(CountingTree owner, ISuffixTreeVisitor inner) : ISuffixTreeVisitor
    {
        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            owner.Hit("ISuffixTreeVisitor.VisitNode(Int32,Int32,Int32,Int32,Int32)");
            inner.VisitNode(startIndex, endIndex, leafCount, childCount, depth);
        }

        public void EnterBranch(int key)
        {
            owner.Hit("ISuffixTreeVisitor.EnterBranch(Int32)");
            inner.EnterBranch(key);
        }

        public void ExitBranch()
        {
            owner.Hit("ISuffixTreeVisitor.ExitBranch()");
            inner.ExitBranch();
        }
    }
}
