using SuffixTree.Indexing;

namespace SuffixTree.Esa;

/// <summary>
/// The enhanced suffix array behind the index contract (work package 1b.5, ADR-001 section 6): <see cref="ISuffixIndex"/>, <see cref="ISuffixTreeNavigation"/> and <see cref="ISuffixArrayAccess"/>
/// over an in-memory <see cref="EsaIndex"/>. The tree is the terminated suffix tree of ADR-001 section 4; a <see cref="IndexLocus"/> is <c>(Lo, Hi, Length, Aux)</c> with <c>[Lo..Hi]</c> the ranks of the
/// suffixes that start with the locus string, <c>Length</c> the number of matched symbols (a leaf counts its terminator) and <c>Aux</c> the depth of the lowest explicit node at or below the locus (so a locus is
/// an explicit node exactly when <c>Length == Aux</c>). The root is <c>(0, n, 0, 0)</c>; the leaf of the empty suffix, the terminator child of the root, is <c>(0, 0, 1, 1)</c> and does not exist for the empty text.
/// <para>
/// Counts and locates ignore rank 0 (the empty suffix is no text position), so <c>Count(Root) == n</c> and the children of a node always add up to its count. Everything that the core index has already -
/// search, count, locate, children, edge symbols, suffix-array and LCP access - runs without the navigation pack; the operations that need ISA or the LCP pyramid (<c>TryParent</c>, <c>TrySuffixLink</c>,
/// <c>TryExtendLeft</c>, <c>RankOf</c>, <c>Lcp</c>, <c>Lca</c>) build it on first use (once, thread-safe, <see cref="EsaIndex.Navigation"/>) and are <see cref="IndexSupport.Degraded"/> until then.
/// Loci are checked against the index: a tuple that cannot be a locus of this text is an <see cref="ArgumentException"/> (cheap range checks everywhere; the navigation operations also read the LCP array
/// at both ends of the interval). Immutable after construction and safe for concurrent readers.
/// </para>
/// </summary>
internal sealed partial class EsaSuffixIndex : ISuffixTreeNavigation, ISuffixArrayAccess
{
    private const string PackCost = "builds ISA and the LCP pyramid on first use, O(n), 4.13 B/char, about 0.1 s at 10 Mbp";
    private const string ExtendLeftCost = "O(width of the locus interval): reads SA, the text and ISA over the interval, then one range minimum; " + PackCost;

    private readonly EsaIndex _index;
    private readonly long _n;
    private readonly IIndexText _text;

    /// <summary>The contract view of <paramref name="index"/> (shared, not copied).</summary>
    /// <exception cref="ArgumentNullException"><paramref name="index"/> is null.</exception>
    internal EsaSuffixIndex(EsaIndex index)
    {
        ArgumentNullException.ThrowIfNull(index);
        _index = index;
        _n = index.TextLength;
        _text = index.TextView.AsIndexText();
        _isMapped = index.IsMapped;
    }

    /// <summary>Builds the index of <paramref name="text"/> (SA-IS, LCP, child table) and wraps it.</summary>
    /// <exception cref="ArgumentNullException"><paramref name="text"/> is null.</exception>
    internal static EsaSuffixIndex Build(string text, ChildLayoutKind? layout = null) => new(EsaIndex.Build(text, layout));

    /// <inheritdoc />
    public string Engine => "esa";

    /// <inheritdoc />
    public long TextLength => _n;

    /// <inheritdoc />
    public IIndexText Text => _text;

    /// <inheritdoc />
    public IndexLocus Root => new(0, _n, 0, 0);

    /// <summary>
    /// ADR-001 section 6, column ESA: native for count, extend right, locate, text, child, string depth, suffix-array and LCP access; <see cref="IndexOperation.ExtendLeft"/> is always degraded (O(width));
    /// parent, suffix link, ISA, <c>Lcp(i, j)</c> and <c>Lca</c> are degraded until the navigation pack exists and native afterwards.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="operation"/> is not an <see cref="IndexOperation"/>.</exception>
    public IndexSupport Support(IndexOperation operation) => operation switch
    {
        IndexOperation.Count or IndexOperation.ExtendRight or IndexOperation.LazyLocate or IndexOperation.TextAccess
            or IndexOperation.Child or IndexOperation.StringDepth or IndexOperation.SuffixArrayAccess or IndexOperation.LcpAccess => IndexSupport.Native,
        IndexOperation.ExtendLeft => IndexSupport.Degraded,
        IndexOperation.Parent or IndexOperation.SuffixLink or IndexOperation.InverseSuffixArrayAccess or IndexOperation.Lcp or IndexOperation.Lca => _index.SuffixLinkSupport,
        _ => throw new ArgumentOutOfRangeException(nameof(operation), operation, "not an IndexOperation"),
    };

    /// <inheritdoc />
    public string? SupportNote(IndexOperation operation)
    {
        IndexSupport support = Support(operation);
        if (operation == IndexOperation.ExtendLeft) return ExtendLeftCost;
        return support == IndexSupport.Degraded ? PackCost : null;
    }

    /// <inheritdoc />
    public long Count(in IndexLocus locus)
    {
        CheckLocus(locus, nameof(locus));
        return locus.Hi - Math.Max(locus.Lo, 1) + 1;
    }

    /// <inheritdoc />
    public bool TryExtendRight(in IndexLocus locus, int symbol, out IndexLocus result)
    {
        if (symbol is < 0 or > char.MaxValue)
        {
            Resolve(locus, nameof(locus));
            result = default;
            return false;
        }
        ReadOnlySpan<char> one = [(char)symbol];
        return Extend(locus, one, out result);
    }

    /// <inheritdoc />
    public bool TryExtendRight(in IndexLocus locus, ReadOnlySpan<char> pattern, out IndexLocus result) => Extend(locus, pattern, out result);

    /// <inheritdoc />
    public bool TryExtendLeft(in IndexLocus locus, int symbol, out IndexLocus result)
    {
        IndexLocus at = locus; // a copy: result may be the same variable as locus (the contract), and result is written before the locus is read again
        Resolve(at, nameof(locus));
        result = default;
        if (symbol is < 0 or > char.MaxValue) return false;
        if (at.Length == 0) return TryExtendRight(at, symbol, out result); // the empty string: a symbol on the left is a symbol on the right
        EsaNavigationPack pack = _index.Navigation; // the pack is 32-bit until S2.5: it refuses a text above 2^31 - 2 symbols by name
        SaView sa = _index.Sa;
        EsaText text = _index.TextView;
        long low = long.MaxValue, high = -1;
        for (long rank = at.Lo; rank <= at.Hi; rank++)
        {
            // the suffixes that start with symbol + locus string are the predecessors of the occurrences that are preceded by the symbol: a contiguous rank interval
            long position = sa[rank];
            if (position > 0 && text[position - 1] == symbol)
            {
                long target = pack.Isa(IntBridge.Length(position - 1, "EsaSuffixIndex.TryExtendLeft (the inverse suffix array of the navigation pack is 32-bit until S2.5)"));
                if (target < low) low = target;
                if (target > high) high = target;
            }
        }
        if (high < 0) return false;
        long aux = low == high ? _n - sa[low] + 1 : pack.Pyramid.Min(low + 1, high);
        result = new IndexLocus(low, high, at.Length + 1, aux);
        return true;
    }

    /// <inheritdoc />
    public int Locate(in IndexLocus locus, long skip, Span<long> positions)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(skip);
        CheckLocus(locus, nameof(locus));
        long from = Math.Max(locus.Lo, 1);
        if (from > locus.Hi) return 0;
        return _index.Locate(from, locus.Hi, skip, positions);
    }

    /// <summary>
    /// The extension by a pattern: inside an edge the next symbol is the text at the first suffix of the interval; at a node the child whose edge starts with the symbol is found by walking the l-indices.
    /// <paramref name="result"/> may be the same variable as <paramref name="locus"/> (the contract): everything that is needed from the locus is read before the first write to the result.
    /// </summary>
    private bool Extend(in IndexLocus locus, ReadOnlySpan<char> pattern, out IndexLocus result)
    {
        EsaNode64 node = Resolve(locus, nameof(locus));
        long length = locus.Length;
        result = default;
        EsaText text = _index.TextView;
        foreach (char symbol in pattern)
        {
            if (length < node.Depth)
            {
                long at = node.Pos + length;
                if ((ulong)at >= (ulong)_n || text[at] != symbol) return false;
            }
            else
            {
                if (!FindChild(node, symbol, out EsaNode64 next)) return false;
                node = next;
            }
            length++;
        }
        result = new IndexLocus(node.Lb, node.Rb, length, node.Depth);
        return true;
    }

    /// <summary>The O(1) checks that every operation applies: the tuple lies inside the index and has the shape of a locus (root, the leaf of the empty suffix, or a point of a non-empty suffix).</summary>
    private void CheckLocus(in IndexLocus locus, string parameterName)
    {
        long n = _n;
        bool valid = locus.Lo >= 0 && locus.Lo <= locus.Hi && locus.Hi <= n && locus.Length >= 0 && locus.Aux >= locus.Length && locus.Aux <= n + 1;
        if (valid)
        {
            valid = locus.Length == 0
                ? locus.Lo == 0 && locus.Hi == n && locus.Aux == 0
                : n > 0 && (locus.Lo > 0 || (locus.Hi == 0 && locus.Aux == 1));
        }
        if (!valid) throw new ArgumentException($"{locus} is not a locus of this index (text length {n})", parameterName);
    }

    /// <summary>
    /// The node below (or at) a locus, as the internal handle: the interval, its depth (<c>Aux</c>), the depth of its parent and its last-child flag from the LCP values just outside the interval
    /// (<c>parent depth = max(LCP[lb], LCP[rb + 1])</c>; the interval is the last child exactly when <c>LCP[lb] &gt; LCP[rb + 1]</c>), and <c>Pos = SA[lb]</c>. Checks that the interval really is the interval of
    /// that node (a leaf has the depth of its suffix plus the terminator; an inner interval has its first l-index in the child table with the value <c>Aux</c> and smaller LCP values at both ends) and that
    /// the matched length lies on the edge above it.
    /// </summary>
    private EsaNode64 Resolve(in IndexLocus locus, string parameterName)
    {
        CheckLocus(locus, parameterName);
        long depth = locus.Aux;
        if (depth == 0) return new EsaNode64(0, _n, 0, 0, _n, true);
        long lb = locus.Lo;
        long rb = locus.Hi;
        CompactLcp lcp = _index.Lcp;
        long left = lcp[lb];
        long right = rb < _n ? lcp[rb + 1] : -1;
        long parentDepth = Math.Max(left, right);
        bool lastChild = left > right;
        long pos = _index.Sa[lb];
        bool valid = parentDepth < locus.Length;
        if (valid)
        {
            if (lb == rb) valid = depth == _n - pos + 1;
            else
            {
                long first = _index.Child.FirstLIndex(lb, rb, lastChild);
                valid = first > lb && first <= rb && lcp[first] == depth;
            }
        }
        if (!valid) throw new ArgumentException($"{locus} is not a locus of this index (no node of depth {depth} spans ranks {lb}..{rb})", parameterName);
        return new EsaNode64(lb, rb, depth, parentDepth, pos, lastChild);
    }

    /// <summary><see cref="Resolve"/> for the operations that take an explicit node.</summary>
    private EsaNode64 ResolveNode(in IndexLocus node, string parameterName)
    {
        EsaNode64 resolved = Resolve(node, parameterName);
        if (node.Length != node.Aux) throw new ArgumentException($"{node} is a point inside an edge, not an explicit node", parameterName);
        return resolved;
    }

    private static IndexLocus ToLocus(in EsaNode64 node) => new(node.Lb, node.Rb, node.Depth, node.Depth);

    private void CheckRank(long rank, string parameterName)
    {
        if ((ulong)rank > (ulong)_n) throw new ArgumentOutOfRangeException(parameterName, rank, $"not a rank or position of a text of {_n} symbols (0..{_n})");
    }
}
