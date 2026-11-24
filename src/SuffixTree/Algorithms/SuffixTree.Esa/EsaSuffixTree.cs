using System.Runtime.CompilerServices;
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
/// Lifetime: <see cref="Dispose"/> is idempotent; every member of a disposed tree throws <see cref="ObjectDisposedException"/> (before any argument check). Reads are safe for concurrent readers, and so is
/// a <see cref="Dispose"/> that arrives while another thread is inside a call. For a tree over a container (work package 1b.8) every public member registers itself as a reader for the length of the call: since work package 1b.9g
/// by a plain store of the tree's id into the cell of the calling thread (<see cref="ReaderCell"/>, no atomic operation and no fence on the path of a call), and <see cref="Dispose"/> raises a flag, runs a process-wide barrier
/// and reads the cells (<see cref="ReaderGate{TSeams}"/>; the proof is in ADR-004 section 5a, "1b.9g"). The pins of the navigation pack and the mapped container are released by whoever finds the tree disposed and no reader inside, so a call that was in flight when <see cref="Dispose"/> ran
/// completes with its answer, and a call that starts afterwards throws <see cref="ObjectDisposedException"/>; nothing reads the mapping after it was unmapped. For an in-memory tree (work package 1b.9b) there is nothing to wait for:
/// <see cref="Dispose"/> only raises the flag and releases nothing that a reader uses (the arrays are managed memory, the pins of the navigation pack are released by the finalizer of the pack's pinned array when the tree is
/// unreachable), so its members check the flag and take no count, and the same two statements hold. <see cref="Dispose"/> does not wait for the readers (a
/// <see cref="Dispose"/> from inside a visitor of <see cref="Traverse"/> cannot deadlock). The results that read the mapping later (<see cref="Text"/>, <c>ISuffixIndex.Text</c>, the sequence of
/// <see cref="EnumerateSuffixes"/>) of a file-backed tree count themselves as readers at every access. The
/// results the contract leaves open are those of the other backends' contract: the order of <see cref="FindAllOccurrences(string)"/> is suffix-array order, the representative of
/// <see cref="LongestRepeatedSubstring"/> on a tie is the lexicographically smallest tied substring, the text position of <see cref="LongestCommonSubstringInfo(string)"/> and of the anchors is the first
/// occurrence the navigator offers.
/// </remarks>
public sealed partial class EsaSuffixTree : ISuffixTree, ISuffixTreeNavigation, ISuffixArrayAccess, IDisposable
{
    private readonly EsaIndex _index;
    private readonly EsaSuffixIndex _contract;
    private readonly ITextSource _text;
    private readonly EsaStatistics _statistics;
    private int _disposed;
    private readonly ReaderGate<ReaderGateSeams>? _gate;
    private string? _longestRepeat;
    private readonly External.SesaFile? _file;
    private readonly ITextSource? _guardedText;
    private readonly IIndexText? _guardedIndexText;

    /// <summary>Wraps an index. <paramref name="statistics"/> replaces the header statistics of the index (the tests use it to feed a synthetic node count to the guard).</summary>
    /// <exception cref="ArgumentNullException"><paramref name="index"/> is null.</exception>
    /// <exception cref="NotSupportedException">The text is beyond the 32-bit surface (<see cref="EnsureTextLength"/>).</exception>
    internal EsaSuffixTree(EsaIndex index, EsaStatistics? statistics = null)
    {
        ArgumentNullException.ThrowIfNull(index);
        EnsureTextLength(index.TextLength);
        _index = index;
        _statistics = statistics ?? index.Statistics;
        _text = index.TextView.AsTextSource();
        _contract = new EsaSuffixIndex(index);
        _file = index.File;
        if (_file is not null)
        {
            _guardedText = new GuardedTextSource(this, _text);
            _guardedIndexText = new GuardedIndexText(this, _contract.Text);
            _gate = new ReaderGate<ReaderGateSeams>(this, ReleaseResources);
        }
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

    /// <summary>
    /// The tree over an opened container, which it owns: <see cref="Dispose"/> releases the mapping. When the tree cannot be made (a text beyond the 32-bit surface) the file is disposed here and the exception propagates.
    /// </summary>
    /// <exception cref="NotSupportedException">The text is beyond the 32-bit surface (<see cref="EnsureTextLength"/>).</exception>
    internal static EsaSuffixTree FromFile(External.SesaFile file)
    {
        ArgumentNullException.ThrowIfNull(file);
        try
        {
            EnsureTextLength(file.Header.TextLength);
            return new EsaSuffixTree(EsaIndex.FromFile(file));
        }
        catch
        {
            file.Dispose();
            throw;
        }
    }

    /// <summary>The container this tree reads in place, or null for a tree built in memory (1b.7f).</summary>
    internal External.SesaFile? File => _file;

    /// <summary>The id that the calls of this tree carry in the reader cells of their threads (work package 1b.9g), or 0 for a tree built in memory, whose calls are not registered.</summary>
    internal long ReaderId => _gate?.Id ?? 0;

    // ------------------------------------------------------------------ lifetime

    /// <summary>
    /// Ends the life of the tree: new calls throw <see cref="ObjectDisposedException"/>, the calls in flight complete, and for a tree over a container the pins of the navigation pack (if it was built; the disposal never builds
    /// it) and the mapped container are released when the last of them has returned; an in-memory tree releases nothing here (its pins go with the finalizer of the pack's pinned array). Idempotent; does not wait.
    /// </summary>
    public void Dispose()
    {
        if (_gate is not null)
        {
            _gate.Dispose();
            return;
        }

        // an in-memory tree releases nothing here: what a reader of it uses (the arrays, and the pins of the navigation pack, which go with the finalizer of the pack's pinned array) stays valid until the tree is unreachable
        Interlocked.CompareExchange(ref _disposed, 1, 0);
    }

    /// <summary>
    /// Registers a call as a reader, or throws <see cref="ObjectDisposedException"/> when the tree is disposed. A tree over a container goes through its <see cref="ReaderGate{TSeams}"/>: the id of the tree is stored in the cell of the
    /// calling thread and the flag is read after it, and <see cref="Dispose"/> raises the flag, runs a process-wide barrier and reads the cells, so a call is either seen by the releaser or sees the flag and goes away before it has touched
    /// anything (ADR-004 section 5a, "1b.9g"); the path of a call has no atomic operation and no fence. An in-memory tree takes no entry: <see cref="Dispose"/> releases nothing that a call in flight uses (work package 1b.9b), so the call needs the
    /// flag and nothing else; the scope still holds the tree, which keeps it reachable for the length of the call.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private ReadScope Enter()
    {
        ReaderGate<ReaderGateSeams>? gate = _gate;
        if (gate is not null) return new ReadScope(gate, gate.Enter());
        ObjectDisposedException.ThrowIf(Volatile.Read(ref _disposed) != 0, this);
        return new ReadScope(this, null);
    }

    /// <summary>What the last reader out (or <see cref="Dispose"/>) runs once, for a tree over a container: the pins of the navigation pack and the mapping.</summary>
    private void ReleaseResources()
    {
        if (_index.NavigationBuilt) _index.Navigation.Dispose();
        _file?.Dispose();
    }

    /// <summary>
    /// The registration of one call (a using declaration at the top of every public member). It holds the object that keeps the tree reachable for the length of the call: the gate of a tree over a container (which holds the tree)
    /// together with the cell of the thread, or the in-memory tree itself (and no cell): the collector must not finalize the pins of the navigation pack under a reader that reads through pointers while the compiler already
    /// regards <c>this</c> as dead.
    /// </summary>
    private readonly struct ReadScope(object holder, ReaderCell? cell) : IDisposable
    {
        public void Dispose()
        {
            if (cell is not null) Unsafe.As<ReaderGate<ReaderGateSeams>>(holder).Leave(cell);
        }
    }

    // ------------------------------------------------------------------ ISuffixTreeSearch

    /// <inheritdoc />
    public ITextSource Text
    {
        get
        {
            using ReadScope scope = Enter();
            return _guardedText ?? _text;
        }
    }

    /// <inheritdoc />
    /// <exception cref="NotSupportedException">The node count does not fit an <see cref="int"/> (<see cref="CheckedNodeCount"/>).</exception>
    public int NodeCount
    {
        get
        {
            using ReadScope scope = Enter();
            return CheckedNodeCount(_statistics.NodeCount);
        }
    }

    /// <inheritdoc />
    public int LeafCount
    {
        get
        {
            using ReadScope scope = Enter();
            return IntBridge.Length(_statistics.LeafCount, "EsaSuffixTree.LeafCount (the 64-bit member is ISuffixIndex.TextLength)");
        }
    }

    /// <inheritdoc />
    public int MaxDepth
    {
        get
        {
            using ReadScope scope = Enter();
            return IntBridge.Length(_statistics.MaxDepth, "EsaSuffixTree.MaxDepth (the 64-bit member is ISuffixIndex.TextLength)");
        }
    }

    /// <inheritdoc />
    public bool IsEmpty
    {
        get
        {
            using ReadScope scope = Enter();
            return _index.Length == 0;
        }
    }

    /// <inheritdoc />
    public bool Contains(string value)
    {
        using ReadScope scope = Enter();
        SuffixTreeSearchContracts.EnsureNotNull(value, nameof(value));
        return _index.Contains(value.AsSpan());
    }

    /// <inheritdoc />
    public bool Contains(ReadOnlySpan<char> value)
    {
        using ReadScope scope = Enter();
        return _index.Contains(value);
    }

    /// <inheritdoc />
    public IReadOnlyList<int> FindAllOccurrences(string pattern)
    {
        using ReadScope scope = Enter();
        SuffixTreeSearchContracts.EnsureNotNull(pattern, nameof(pattern));
        return _index.FindAllOccurrences(pattern.AsSpan());
    }

    /// <inheritdoc />
    public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern)
    {
        using ReadScope scope = Enter();
        return _index.FindAllOccurrences(pattern);
    }

    /// <inheritdoc />
    public int CountOccurrences(string pattern)
    {
        using ReadScope scope = Enter();
        SuffixTreeSearchContracts.EnsureNotNull(pattern, nameof(pattern));
        return _index.CountOccurrences(pattern.AsSpan());
    }

    /// <inheritdoc />
    public int CountOccurrences(ReadOnlySpan<char> pattern)
    {
        using ReadScope scope = Enter();
        return _index.CountOccurrences(pattern);
    }

    // ------------------------------------------------------------------ ISuffixTreeAnalysis and the members of ISuffixTree

    /// <inheritdoc />
    /// <remarks>The string is computed once and kept, as the other backends do (repeated calls return the same instance); concurrent first calls may both compute it, the result is the same string.</remarks>
    public string LongestRepeatedSubstring()
    {
        using ReadScope scope = Enter();
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
        using ReadScope scope = Enter();
        return _index.LongestRepeatedSubstringMemory();
    }

    /// <inheritdoc />
    public IReadOnlyList<string> GetAllSuffixes()
    {
        using ReadScope scope = Enter();
        return _index.GetAllSuffixes();
    }

    /// <inheritdoc />
    /// <remarks>The disposal check is eager, the sequence itself is lazy (as for the persistent tree).</remarks>
    public IEnumerable<string> EnumerateSuffixes()
    {
        using ReadScope scope = Enter();
        IEnumerable<string> sequence = _index.EnumerateSuffixes();
        return _file is null ? sequence : new GuardedSequence<string>(this, sequence);
    }

    /// <inheritdoc />
    public string LongestCommonSubstring(string other)
    {
        using ReadScope scope = Enter();
        return _index.LongestCommonSubstring(other);
    }

    /// <inheritdoc />
    public string LongestCommonSubstring(ReadOnlySpan<char> other)
    {
        using ReadScope scope = Enter();
        return _index.LongestCommonSubstring(other);
    }

    /// <inheritdoc />
    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other)
    {
        using ReadScope scope = Enter();
        return _index.LongestCommonSubstringInfo(other);
    }

    /// <summary>The same for a character span (a member of the class, not of the interface, as for the persistent tree).</summary>
    /// <param name="other">The other string.</param>
    /// <returns>The longest common substring with one position in the text and its first position in <paramref name="other"/>; <c>("", -1, -1)</c> when there is none.</returns>
    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(ReadOnlySpan<char> other)
    {
        using ReadScope scope = Enter();
        return _index.LongestCommonSubstringInfo(other);
    }

    /// <inheritdoc />
    public (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other)
    {
        using ReadScope scope = Enter();
        return _index.FindAllLongestCommonSubstrings(other);
    }

    /// <inheritdoc />
    public IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)> FindAllDistinctLongestCommonSubstrings(string other)
    {
        using ReadScope scope = Enter();
        return _index.FindAllDistinctLongestCommonSubstrings(other);
    }

    /// <inheritdoc />
    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(string query, int minLength)
    {
        using ReadScope scope = Enter();
        return _index.FindExactMatchAnchors(query, minLength);
    }

    /// <inheritdoc />
    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindMaximalExactMatches(string query, int minLength)
    {
        using ReadScope scope = Enter();
        return _index.FindMaximalExactMatches(query, minLength);
    }

    /// <inheritdoc />
    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindMaximalUniqueMatches(string query, int minLength, MumUniqueness uniqueness = MumUniqueness.Both)
    {
        using ReadScope scope = Enter();
        return _index.FindMaximalUniqueMatches(query, minLength, uniqueness);
    }

    /// <inheritdoc />
    public long[] CountDistinctSubstringsByLength(int maxLength)
    {
        using ReadScope scope = Enter();
        return _index.CountDistinctSubstringsByLength(maxLength);
    }

    /// <inheritdoc />
    public long CountDistinctSubstrings()
    {
        using ReadScope scope = Enter();
        return _index.CountDistinctSubstrings();
    }

    /// <inheritdoc />
    public IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> FindAllLongestRepeatedSubstrings()
    {
        using ReadScope scope = Enter();
        return _index.FindAllLongestRepeatedSubstrings();
    }

    /// <inheritdoc />
    public IReadOnlyList<(int FirstPosition, int SecondPosition, int Length)> FindMaximalRepeatedPairs(int minLength, Func<char, bool>? isUniqueSymbol = null)
    {
        using ReadScope scope = Enter();
        return SuffixTreeAlgorithms.FindMaximalRepeatedPairs(this, minLength, isUniqueSymbol);
    }

    // ------------------------------------------------------------------ results that read the mapping after the call has returned

    /// <summary>The text of a file-backed tree: every access is a call of the tree and counts as a reader (a stored copy of the text must not read the mapping after <see cref="Dispose"/>).</summary>
    private sealed class GuardedTextSource(EsaSuffixTree owner, ITextSource inner) : ITextSource
    {
        public int Length => inner.Length;

        public char this[int index]
        {
            get
            {
                using ReadScope scope = owner.Enter();
                return inner[index];
            }
        }

        public string Substring(int start, int length)
        {
            using ReadScope scope = owner.Enter();
            return inner.Substring(start, length);
        }

        public ReadOnlySpan<char> Slice(int start, int length)
        {
            using ReadScope scope = owner.Enter();
            return inner.Slice(start, length);
        }

        public IEnumerator<char> GetEnumerator() => new GuardedEnumerator<char>(owner, inner.GetEnumerator());

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() => GetEnumerator();

        public override string ToString()
        {
            using ReadScope scope = owner.Enter();
            return inner.ToString() ?? string.Empty;
        }
    }

    /// <summary>The <see cref="IIndexText"/> of a file-backed tree, guarded as <see cref="GuardedTextSource"/>.</summary>
    private sealed class GuardedIndexText(EsaSuffixTree owner, IIndexText inner) : IIndexText
    {
        public long Length => inner.Length;

        public int SymbolAt(long position)
        {
            using ReadScope scope = owner.Enter();
            return inner.SymbolAt(position);
        }

        public int CopyTo(long start, Span<char> destination)
        {
            using ReadScope scope = owner.Enter();
            return inner.CopyTo(start, destination);
        }
    }

    private sealed class GuardedSequence<T>(EsaSuffixTree owner, IEnumerable<T> inner) : IEnumerable<T>
    {
        public IEnumerator<T> GetEnumerator() => new GuardedEnumerator<T>(owner, inner.GetEnumerator());

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() => GetEnumerator();
    }

    /// <summary>An enumerator whose every step is a call of the tree (the lazy walk reads the mapping at each <c>MoveNext</c>).</summary>
    private sealed class GuardedEnumerator<T>(EsaSuffixTree owner, IEnumerator<T> inner) : IEnumerator<T>
    {
        public T Current => inner.Current;

        object? System.Collections.IEnumerator.Current => Current;

        public bool MoveNext()
        {
            using ReadScope scope = owner.Enter();
            return inner.MoveNext();
        }

        public void Reset() => inner.Reset();

        public void Dispose() => inner.Dispose();
    }
}
