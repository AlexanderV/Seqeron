using SuffixTree.Indexing;

namespace SuffixTree.Esa;

/// <summary>
/// The enhanced suffix array as an <see cref="ISuffixTree"/> (work package 1b.6, ADR-001 section 7): a facade with the <see cref="int"/> surface of the existing trees over the in-memory
/// <see cref="EsaIndex"/>, and with the <c>long</c> navigation and suffix-array surfaces of the index contract (<see cref="ISuffixTreeNavigation"/>, <see cref="ISuffixArrayAccess"/>) over the same
/// object. The tree is node for node the terminated suffix tree of the text: node counts, the <see cref="ISuffixTreeDiagnostics.Traverse"/> event stream and the logical hash equal those of the
/// Ukkonen tree and of the persistent v6 tree. Search, the native analysis operations (longest repeat, distinct substrings, suffixes) and <c>Traverse</c> need only the suffix array, the LCP array and the
/// child table; the operations that need suffix links (the longest common substring family, anchors, maximal exact and unique matches, <see cref="PrintTree"/>) and the navigation members that
/// need the inverse suffix array build the navigation pack once, on first use (<see cref="IndexSupport"/> reports <see cref="IndexSupport.Degraded"/> until then).
/// </summary>
/// <remarks>
/// Lifetime is that of <c>PersistentSuffixTree</c>: <see cref="Dispose"/> is idempotent and releases the pins of the navigation pack; every member of a disposed tree throws
/// <see cref="ObjectDisposedException"/> (before any argument check). Reads are safe for concurrent readers; <see cref="Dispose"/> must not run while another thread is inside a call. The
/// results the contract leaves open are those of the other backends' contract: the order of <see cref="FindAllOccurrences(string)"/> is suffix-array order, the representative of
/// <see cref="LongestRepeatedSubstring"/> on a tie is the lexicographically smallest tied substring, the text position of <see cref="LongestCommonSubstringInfo(string)"/> and of the anchors is the first
/// occurrence the navigator offers.
/// </remarks>
public sealed partial class EsaSuffixTree : ISuffixTree, ISuffixTreeNavigation, ISuffixArrayAccess, IDisposable
{
    private readonly EsaIndex _index;
    private readonly EsaSuffixIndex _contract;
    private readonly StringTextSource _text;
    private readonly EsaStatistics _statistics;
    private int _disposed;
    private string? _longestRepeat;

    /// <summary>Wraps an index. <paramref name="statistics"/> replaces the header statistics of the index (the tests use it to feed a synthetic node count to the guard).</summary>
    /// <exception cref="ArgumentNullException"><paramref name="index"/> is null.</exception>
    /// <exception cref="NotSupportedException">The text is beyond the 32-bit surface (<see cref="EnsureTextLength"/>).</exception>
    internal EsaSuffixTree(EsaIndex index, EsaStatistics? statistics = null)
    {
        ArgumentNullException.ThrowIfNull(index);
        EnsureTextLength(index.Length);
        _index = index;
        _statistics = statistics ?? index.Statistics;
        _text = new StringTextSource(index.Text);
        _contract = new EsaSuffixIndex(index);
    }

    /// <summary>Builds the tree of <paramref name="text"/> with the default child-table layout.</summary>
    /// <param name="text">The text; the tree shares the string, it is not copied.</param>
    /// <returns>The tree.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="text"/> is null.</exception>
    public static EsaSuffixTree Build(string text) => Build(text, null);

    internal static EsaSuffixTree Build(string text, ChildLayoutKind? layout)
    {
        ArgumentNullException.ThrowIfNull(text);
        return new EsaSuffixTree(EsaIndex.Build(text, layout));
    }

    /// <summary>The largest text the 32-bit <see cref="ISuffixTree"/> surface carries: the root of <see cref="ISuffixTreeDiagnostics.Traverse"/> reports <c>n + 1</c> leaves as an <see cref="int"/>.</summary>
    internal const long MaxTextLength = int.MaxValue - 1L;

    /// <summary>
    /// The guard of ADR-001 sections 7 and 8 for a text of <paramref name="textLength"/> symbols (a container-opened text can be longer than a managed string): a text beyond
    /// <see cref="MaxTextLength"/> is refused with a clear exception instead of an overflow.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">The length is negative.</exception>
    /// <exception cref="NotSupportedException">The length is above <see cref="MaxTextLength"/>.</exception>
    internal static void EnsureTextLength(long textLength)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(textLength);
        if (textLength > MaxTextLength)
            throw new NotSupportedException($"A text of {textLength} symbols is beyond the 32-bit ISuffixTree surface (at most {MaxTextLength}); use the long API (ISuffixIndex).");
    }

    /// <summary>The node count as an <see cref="int"/>, or a clear exception when it does not fit: never a wrapped value.</summary>
    /// <exception cref="NotSupportedException"><paramref name="nodeCount"/> is above <see cref="int.MaxValue"/>.</exception>
    internal static int CheckedNodeCount(long nodeCount)
    {
        if (nodeCount > int.MaxValue)
            throw new NotSupportedException($"The suffix tree has {nodeCount} nodes, more than the 32-bit ISuffixTree surface can report ({int.MaxValue}); use the 64-bit ISuffixIndex surface.");
        return (int)nodeCount;
    }

    /// <summary>The wrapped index (tests inspect the navigation pack through it); not guarded by the disposal check.</summary>
    internal EsaIndex Index => _index;

    // ------------------------------------------------------------------ lifetime

    /// <summary>Releases the pins of the navigation pack, if it was built (the disposal never builds it). Idempotent.</summary>
    public void Dispose()
    {
        if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0) return;
        if (_index.NavigationBuilt) _index.Navigation.Dispose();
    }

    private void ThrowIfDisposed() => ObjectDisposedException.ThrowIf(Volatile.Read(ref _disposed) != 0, this);

    // ------------------------------------------------------------------ ISuffixTreeSearch

    /// <inheritdoc />
    public ITextSource Text
    {
        get
        {
            ThrowIfDisposed();
            return _text;
        }
    }

    /// <inheritdoc />
    /// <exception cref="NotSupportedException">The node count does not fit an <see cref="int"/> (<see cref="CheckedNodeCount"/>).</exception>
    public int NodeCount
    {
        get
        {
            ThrowIfDisposed();
            return CheckedNodeCount(_statistics.NodeCount);
        }
    }

    /// <inheritdoc />
    public int LeafCount
    {
        get
        {
            ThrowIfDisposed();
            return _statistics.LeafCount;
        }
    }

    /// <inheritdoc />
    public int MaxDepth
    {
        get
        {
            ThrowIfDisposed();
            return _statistics.MaxDepth;
        }
    }

    /// <inheritdoc />
    public bool IsEmpty
    {
        get
        {
            ThrowIfDisposed();
            return _index.Length == 0;
        }
    }

    /// <inheritdoc />
    public bool Contains(string value)
    {
        ThrowIfDisposed();
        SuffixTreeSearchContracts.EnsureNotNull(value, nameof(value));
        return _index.Contains(value.AsSpan());
    }

    /// <inheritdoc />
    public bool Contains(ReadOnlySpan<char> value)
    {
        ThrowIfDisposed();
        return _index.Contains(value);
    }

    /// <inheritdoc />
    public IReadOnlyList<int> FindAllOccurrences(string pattern)
    {
        ThrowIfDisposed();
        SuffixTreeSearchContracts.EnsureNotNull(pattern, nameof(pattern));
        return _index.FindAllOccurrences(pattern.AsSpan());
    }

    /// <inheritdoc />
    public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern)
    {
        ThrowIfDisposed();
        return _index.FindAllOccurrences(pattern);
    }

    /// <inheritdoc />
    public int CountOccurrences(string pattern)
    {
        ThrowIfDisposed();
        SuffixTreeSearchContracts.EnsureNotNull(pattern, nameof(pattern));
        return _index.CountOccurrences(pattern.AsSpan());
    }

    /// <inheritdoc />
    public int CountOccurrences(ReadOnlySpan<char> pattern)
    {
        ThrowIfDisposed();
        return _index.CountOccurrences(pattern);
    }

    // ------------------------------------------------------------------ ISuffixTreeAnalysis and the members of ISuffixTree

    /// <inheritdoc />
    /// <remarks>The string is computed once and kept, as the other backends do (repeated calls return the same instance); concurrent first calls may both compute it, the result is the same string.</remarks>
    public string LongestRepeatedSubstring()
    {
        ThrowIfDisposed();
        return Volatile.Read(ref _longestRepeat) ?? PublishLongestRepeat();
    }

    private string PublishLongestRepeat()
    {
        string computed = _index.LongestRepeatedSubstring();
        return Interlocked.CompareExchange(ref _longestRepeat, computed, null) ?? computed;
    }

    /// <inheritdoc />
    public ReadOnlyMemory<char> LongestRepeatedSubstringMemory()
    {
        ThrowIfDisposed();
        return _index.LongestRepeatedSubstringMemory();
    }

    /// <inheritdoc />
    public IReadOnlyList<string> GetAllSuffixes()
    {
        ThrowIfDisposed();
        return _index.GetAllSuffixes();
    }

    /// <inheritdoc />
    /// <remarks>The disposal check is eager, the sequence itself is lazy (as for the persistent tree).</remarks>
    public IEnumerable<string> EnumerateSuffixes()
    {
        ThrowIfDisposed();
        return _index.EnumerateSuffixes();
    }

    /// <inheritdoc />
    public string LongestCommonSubstring(string other)
    {
        ThrowIfDisposed();
        return _index.LongestCommonSubstring(other);
    }

    /// <inheritdoc />
    public string LongestCommonSubstring(ReadOnlySpan<char> other)
    {
        ThrowIfDisposed();
        return _index.LongestCommonSubstring(other);
    }

    /// <inheritdoc />
    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other)
    {
        ThrowIfDisposed();
        return _index.LongestCommonSubstringInfo(other);
    }

    /// <summary>The same for a character span (a member of the class, not of the interface, as for the persistent tree).</summary>
    /// <param name="other">The other string.</param>
    /// <returns>The longest common substring with one position in the text and its first position in <paramref name="other"/>; <c>("", -1, -1)</c> when there is none.</returns>
    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(ReadOnlySpan<char> other)
    {
        ThrowIfDisposed();
        return _index.LongestCommonSubstringInfo(other);
    }

    /// <inheritdoc />
    public (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other)
    {
        ThrowIfDisposed();
        return _index.FindAllLongestCommonSubstrings(other);
    }

    /// <inheritdoc />
    public IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)> FindAllDistinctLongestCommonSubstrings(string other)
    {
        ThrowIfDisposed();
        return _index.FindAllDistinctLongestCommonSubstrings(other);
    }

    /// <inheritdoc />
    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(string query, int minLength)
    {
        ThrowIfDisposed();
        return _index.FindExactMatchAnchors(query, minLength);
    }

    /// <inheritdoc />
    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindMaximalExactMatches(string query, int minLength)
    {
        ThrowIfDisposed();
        return _index.FindMaximalExactMatches(query, minLength);
    }

    /// <inheritdoc />
    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindMaximalUniqueMatches(string query, int minLength, MumUniqueness uniqueness = MumUniqueness.Both)
    {
        ThrowIfDisposed();
        return _index.FindMaximalUniqueMatches(query, minLength, uniqueness);
    }

    /// <inheritdoc />
    public long[] CountDistinctSubstringsByLength(int maxLength)
    {
        ThrowIfDisposed();
        return _index.CountDistinctSubstringsByLength(maxLength);
    }

    /// <inheritdoc />
    public long CountDistinctSubstrings()
    {
        ThrowIfDisposed();
        return _index.CountDistinctSubstrings();
    }

    /// <inheritdoc />
    public IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> FindAllLongestRepeatedSubstrings()
    {
        ThrowIfDisposed();
        return _index.FindAllLongestRepeatedSubstrings();
    }

    /// <inheritdoc />
    public IReadOnlyList<(int FirstPosition, int SecondPosition, int Length)> FindMaximalRepeatedPairs(int minLength, Func<char, bool>? isUniqueSymbol = null)
    {
        ThrowIfDisposed();
        return SuffixTreeAlgorithms.FindMaximalRepeatedPairs(this, minLength, isUniqueSymbol);
    }
}
