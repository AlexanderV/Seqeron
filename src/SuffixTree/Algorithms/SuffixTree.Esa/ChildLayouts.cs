namespace SuffixTree.Esa;

/// <summary>Layout of the child table (prototype P3, ADR-001 section 5).</summary>
internal enum ChildLayoutKind
{
    /// <summary>Three plain int arrays <c>next</c>, <c>up</c>, <c>down</c> (12 bytes per rank): the reference, no LCP comparison on a read.</summary>
    Reference3,

    /// <summary>One int array (4 bytes per rank); the kind of a stored field is told apart by comparing LCP values on every read (Abouelhoda, Kurtz, Ohlebusch 2004).</summary>
    Packed,

    /// <summary>One uint array (4 bytes per rank) with the kind in the two high bits; reads need no LCP comparison; texts up to 2^30 - 1 symbols.</summary>
    Tagged,
}

/// <summary>
/// Read side of a child-table layout. The children of an lcp-interval <c>[lb..rb]</c> with value <c>ell</c> start at <c>lb</c> and at every l-index, an
/// l-index being a rank <c>k</c> in <c>(lb, rb]</c> with <c>LCP[k] = ell</c>; the layout stores the first l-index of every internal interval and, for every
/// l-index, the next one.
/// </summary>
internal interface IChildTable
{
    /// <summary>
    /// The first l-index of the internal interval <c>[lb..rb]</c>. <paramref name="lastChild"/> must say whether the interval is the last child of its parent
    /// (always true for the root): the field sits at <c>rb</c> (the interval is followed by a sibling or a new parent) or at <c>lb</c>.
    /// </summary>
    int FirstLIndex(int lb, int rb, bool lastChild);

    /// <summary>The l-index after <paramref name="k"/> in the interval that ends at <paramref name="rb"/>, or -1 if <paramref name="k"/> is the last one.</summary>
    int NextLIndex(int k, int rb);
}

/// <summary>
/// Three plain arrays, indexed by rank (0 = undefined; every stored value is at least 1): <c>Next[k]</c> the l-index after <c>k</c>; <c>Up[j]</c> the first
/// l-index of the non-last child that ends at <c>j</c>; <c>Down[i]</c> the first l-index of the last child that starts at <c>i</c>. The write side is the sink of the scan.
/// </summary>
internal readonly struct Reference3Child : IChildTable, IIntervalSink
{
    private readonly int[] _next;
    private readonly int[] _up;
    private readonly int[] _down;

    private Reference3Child(int[] next, int[] up, int[] down)
    {
        _next = next;
        _up = up;
        _down = down;
    }

    internal static Reference3Child Create(int n) => new(new int[n + 1], new int[n + 1], new int[n + 1]);

    internal long PayloadBytes => 4L * (_next.Length + _up.Length + _down.Length);

    public void LIndex(int previous, int k)
    {
        if (previous >= 0) _next[previous] = k;
    }

    public void Close(int ell, int lb, int rb, int first, bool lastChild)
    {
        if (lastChild) _down[lb] = first;
        else _up[rb] = first;
    }

    public int FirstLIndex(int lb, int rb, bool lastChild) => lastChild ? _down[lb] : _up[rb];

    public int NextLIndex(int k, int rb)
    {
        int v = _next[k];
        return v == 0 ? -1 : v;
    }
}

/// <summary>
/// One int array, slot <c>k</c> holds at most one field: the next l-index of <c>k</c> when it is not <c>k + 1</c> (implied by <c>LCP[k + 1] = LCP[k]</c>), else the
/// first l-index of the non-last child that ends at <c>k</c> (when <c>LCP[k] &gt; LCP[k + 1]</c>) or of the last child that starts at <c>k</c> (when <c>LCP[k] &lt; LCP[k + 1]</c> and
/// there is no next l-index). The next field and the up field cannot meet (a next field needs <c>LCP[k] &lt;= LCP[k + 1]</c>, an up field <c>LCP[k] &gt; LCP[k + 1]</c>), nor the next and
/// the down field (a last child has no next l-index); a stored value tells next from down by its LCP: next has <c>LCP[v] = LCP[k]</c>, down has <c>LCP[v] &gt; LCP[k]</c>.
/// </summary>
internal readonly struct PackedChild : IChildTable, IIntervalSink
{
    private readonly int[] _cld;
    private readonly CompactLcp _lcp;

    private PackedChild(int[] cld, CompactLcp lcp)
    {
        _cld = cld;
        _lcp = lcp;
    }

    internal static PackedChild Create(int n, CompactLcp lcp) => new(new int[n + 1], lcp);

    internal long PayloadBytes => 4L * _cld.Length;

    public void LIndex(int previous, int k)
    {
        if (previous >= 0 && k != previous + 1) _cld[previous] = k;
    }

    public void Close(int ell, int lb, int rb, int first, bool lastChild) => _cld[lastChild ? lb : rb] = first;

    public int FirstLIndex(int lb, int rb, bool lastChild) => _cld[lastChild ? lb : rb];

    public int NextLIndex(int k, int rb)
    {
        if (k >= rb) return -1;
        int lk = _lcp[k];
        if (_lcp[k + 1] == lk) return k + 1;
        int v = _cld[k];
        return _lcp[v] == lk ? v : -1;
    }
}

/// <summary>
/// One uint array with the kind in bits 30-31 (0 empty, 1 next, 2 up, 3 down) and the value in bits 0-29: a read needs no LCP comparison, the price is the
/// limit <see cref="MaxLength"/> on the text length (the same bound as the Ukkonen tree, <see cref="SuffixTreeLimits.MaxTextLength"/>).
/// </summary>
internal readonly struct TaggedChild : IChildTable, IIntervalSink
{
    /// <summary>The largest text length whose ranks fit the 30 value bits.</summary>
    internal const int MaxLength = (1 << 30) - 1;

    private const uint ValueMask = (1u << 30) - 1;
    private const uint NextTag = 1u << 30;
    private const uint UpTag = 2u << 30;
    private const uint DownTag = 3u << 30;

    private readonly uint[] _cld;

    private TaggedChild(uint[] cld) => _cld = cld;

    /// <exception cref="NotSupportedException">The text is longer than <see cref="MaxLength"/>.</exception>
    internal static TaggedChild Create(int n)
    {
        if (n > MaxLength)
            throw new NotSupportedException($"The tagged child table holds ranks in 30 bits; a text of {n} symbols needs the plain layout.");
        return new TaggedChild(new uint[n + 1]);
    }

    internal long PayloadBytes => 4L * _cld.Length;

    public void LIndex(int previous, int k)
    {
        if (previous >= 0) _cld[previous] = NextTag | (uint)k;
    }

    public void Close(int ell, int lb, int rb, int first, bool lastChild)
    {
        if (lastChild) _cld[lb] = DownTag | (uint)first;
        else _cld[rb] = UpTag | (uint)first;
    }

    public int FirstLIndex(int lb, int rb, bool lastChild) => (int)(_cld[lastChild ? lb : rb] & ValueMask);

    public int NextLIndex(int k, int rb)
    {
        uint c = _cld[k];
        return (c & ~ValueMask) == NextTag ? (int)(c & ValueMask) : -1;
    }
}
