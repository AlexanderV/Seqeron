using System.Runtime.CompilerServices;

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

    internal static uint PackNext(int value) => NextTag | (uint)value;

    internal static uint PackUp(int value) => UpTag | (uint)value;

    internal static uint PackDown(int value) => DownTag | (uint)value;

    /// <summary>The value in bits 0-29 of a stored field.</summary>
    internal static int ValueOf(uint field) => (int)(field & ValueMask);

    /// <summary>True when a stored field is a next-l-index link (tag 1); up, down and empty fields are not.</summary>
    internal static bool IsNextLink(uint field) => (field & ~ValueMask) == NextTag;

    public void LIndex(int previous, int k)
    {
        if (previous >= 0) _cld[previous] = PackNext(k);
    }

    public void Close(int ell, int lb, int rb, int first, bool lastChild)
    {
        if (lastChild) _cld[lb] = PackDown(first);
        else _cld[rb] = PackUp(first);
    }

    public int FirstLIndex(int lb, int rb, bool lastChild) => ValueOf(_cld[lastChild ? lb : rb]);

    public int NextLIndex(int k, int rb)
    {
        uint c = _cld[k];
        return IsNextLink(c) ? ValueOf(c) : -1;
    }
}

/// <summary>
/// The child table of a SESA container read in place (work package 1b.7f; the tagged codec: 1b.9a; widened by S2.3): CHILD codec 1 is the layout of <see cref="PackedChild"/> (one entry per rank, the kind of a stored field told apart by comparing
/// LCP values), codec 2 is the layout of <see cref="TaggedChild"/> (the kind in bits 30-31, no LCP read), over the mapped section instead of an array. The entries are 4 bytes (<c>POS_WIDTH</c> 4) or, for the packed codec, 5 bytes
/// (<c>POS_WIDTH</c> 5, ADR-005 section 6.1: a tagged table does not exist at width 5). The entries are untrusted. With <c>isChecked</c> (the default access mode) a read of a
/// cell outside the section is an <see cref="InvalidDataException"/>; a value that a read returns is only ever used as a rank, and the readers that follow it (<see cref="CompactLcp"/>, the suffix array view) reject a rank outside
/// their arrays; the walks that call <see cref="NextLIndex(int, int)"/> check that the l-indices increase and stay inside the interval, which bounds them. The pattern descent has a loop of its own for each codec (<see cref="MappedSearch"/>);
/// this struct serves the navigator, the traversal and the contract members. The <see cref="int"/> members are what the 32-bit engine of S1 reads (unchanged on a 4-byte section of at most int.MaxValue cells, the checked bridge on a 5-byte one and on a larger one); the <see cref="long"/> members are
/// the 64-bit reads over <see cref="ChildCells{TPos}"/>.
/// </summary>
internal readonly unsafe struct MappedChild : IChildTable
{
    private const int TaggedValueMask = (1 << 30) - 1;
    private const uint TaggedKindMask = 3u << 30;
    private const uint TaggedNext = 1u << 30;

    private readonly int* _cld;
    private readonly long _count;
    private readonly CompactLcp _lcp;
    private readonly bool _checked;
    private readonly bool _tagged;
    private readonly bool _five;

    /// <summary>The cells that the 32-bit members read through the checked bridge: 5-byte cells, and 4-byte cells of a table of more than <see cref="int.MaxValue"/> cells (a cell may be a rank above <see cref="int.MaxValue"/>).</summary>
    private readonly bool _bridged;

    /// <summary>A table of 4-byte cells (the form of S1).</summary>
    internal MappedChild(int* cld, long count, CompactLcp lcp, bool isChecked, bool tagged = false)
        : this(cld, count, lcp, isChecked, tagged, five: false)
    {
    }

    private MappedChild(int* cld, long count, CompactLcp lcp, bool isChecked, bool tagged, bool five)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(count);
        _cld = cld;
        _count = count;
        _lcp = lcp;
        _checked = isChecked;
        _tagged = tagged;
        _five = five;
        _bridged = five || count > int.MaxValue;
    }

    /// <summary>A table of <paramref name="count"/> cells of <paramref name="posWidth"/> (4 or 5) bytes at <paramref name="cells"/>; a 5-byte table has the 3 bytes of padding of its section and is packed.</summary>
    /// <exception cref="ArgumentOutOfRangeException">The width is not 4 or 5, or the count is negative.</exception>
    /// <exception cref="ArgumentException">A tagged table over 5-byte cells.</exception>
    internal static MappedChild FromMappedWide(byte* cells, long count, int posWidth, CompactLcp lcp, bool isChecked, bool tagged)
    {
        ArgumentNullException.ThrowIfNull(lcp);
        if (posWidth is not (4 or 5)) throw new ArgumentOutOfRangeException(nameof(posWidth), posWidth, "the stored width of a cell is 4 or 5 bytes");
        if (tagged && posWidth == 5) throw new ArgumentException("the tagged codec (CHILD codec 2) exists over 4-byte cells only (ADR-005 section 6.1)", nameof(tagged));
        return new MappedChild((int*)cells, count, lcp, isChecked, tagged, posWidth == 5);
    }

    /// <summary>Number of cells (<c>n + 1</c>).</summary>
    internal long Length => _count;

    /// <summary>True for CHILD codec 2 (the kind of a field in bits 30-31).</summary>
    internal bool IsTagged => _tagged;

    /// <summary>True when the cells are 5 bytes, or a table of more than <see cref="int.MaxValue"/> cells: the 32-bit members read them through the checked bridge and the 32-bit navigator does not serve this table.</summary>
    internal bool IsBridged => _bridged;

    /// <summary>The first cell (the pattern descent of <see cref="MappedSearch"/> reads the section through the pointer).</summary>
    /// <exception cref="NotSupportedException">The cells are 5 bytes: there is no pointer to 32-bit cells.</exception>
    internal int* Pointer => _five ? throw new NotSupportedException("the child table has 5-byte cells (POS_WIDTH 5): there is no pointer to 32-bit cells; use the 64-bit members") : _cld;

    /// <summary>Bytes of the section (a 5-byte table includes the 3 bytes of padding).</summary>
    internal long PayloadBytes => _five ? (5 * _count) + PosReader.Padding40 : 4 * _count;

    public int FirstLIndex(int lb, int rb, bool lastChild)
    {
        if (_bridged) return IntBridge.Value(FirstLIndex((long)lb, rb, lastChild), "a first l-index");
        int at = lastChild ? lb : rb;
        if (_checked && (ulong)at >= (ulong)_count) ThrowOutside(at);
        int cell = _cld[at];
        return _tagged ? cell & TaggedValueMask : cell;
    }

    public int NextLIndex(int k, int rb)
    {
        if (_bridged)
        {
            long next = NextLIndex((long)k, rb);
            return next < 0 ? -1 : IntBridge.Value(next, "an l-index");
        }
        if (_tagged)
        {
            if (_checked && (ulong)k >= (ulong)_count) ThrowOutside(k);
            uint cell = (uint)_cld[k];
            return (cell & TaggedKindMask) == TaggedNext ? (int)(cell & TaggedValueMask) : -1;
        }
        if (k >= rb) return -1;
        if (_checked && (ulong)k >= (ulong)_count) ThrowOutside(k);
        int lk = _lcp[k];
        if (_lcp[k + 1] == lk) return k + 1;
        int v = _cld[k];
        return _lcp[v] == lk ? v : -1;
    }

    /// <summary>The first l-index of the interval <c>[lb..rb]</c> at 64-bit ranks (any stored width).</summary>
    /// <exception cref="InvalidDataException">The cell is outside the section (checked mode).</exception>
    internal long FirstLIndex(long lb, long rb, bool lastChild)
    {
        if (_five) return _checked ? Cells(new Pos40<CheckedAccess>((byte*)_cld, _count)).FirstLIndex(lb, rb, lastChild) : Cells(new Pos40<TrustedAccess>((byte*)_cld, _count)).FirstLIndex(lb, rb, lastChild);
        return _checked ? Cells(new Pos32<CheckedAccess>((byte*)_cld, _count)).FirstLIndex(lb, rb, lastChild) : Cells(new Pos32<TrustedAccess>((byte*)_cld, _count)).FirstLIndex(lb, rb, lastChild);
    }

    /// <summary>The l-index after <paramref name="k"/> in the interval that ends at <paramref name="rb"/>, or -1, at 64-bit ranks (any stored width).</summary>
    /// <exception cref="InvalidDataException">A cell or an LCP value is outside its section (checked mode), or the walk meets a rank that is no rank.</exception>
    internal long NextLIndex(long k, long rb)
    {
        if (_five) return _checked ? Cells(new Pos40<CheckedAccess>((byte*)_cld, _count)).NextLIndex(k, rb) : Cells(new Pos40<TrustedAccess>((byte*)_cld, _count)).NextLIndex(k, rb);
        return _checked ? Cells(new Pos32<CheckedAccess>((byte*)_cld, _count)).NextLIndex(k, rb) : Cells(new Pos32<TrustedAccess>((byte*)_cld, _count)).NextLIndex(k, rb);
    }

    private ChildCells<TPos> Cells<TPos>(TPos reader)
        where TPos : struct, IPosReader => new(reader, _lcp, _tagged);

    private void ThrowOutside(long cell) => throw new InvalidDataException($"the child table is read at cell {cell}, but it has {_count} cells: the file is damaged");
}

/// <summary>
/// The child cells over a position reader (ADR-005 section 3.3, work package S2.3): the generic form of <see cref="MappedChild"/> with 64-bit ranks. <c>tagged</c> (CHILD codec 2: the kind of a field in bits 30-31) exists only over a 4-byte
/// reader; the packed cells (codec 1) over any. The checks are the reader's (Checked mode bounds the cell) and the LCP array's (a rank that a cell sent the walk to is checked against the array).
/// </summary>
internal readonly struct ChildCells<TPos>
    where TPos : struct, IPosReader
{
    private const ulong TaggedValueMask = (1UL << 30) - 1;
    private const ulong TaggedKindMask = 3UL << 30;
    private const ulong TaggedNext = 1UL << 30;

    private readonly TPos _reader;
    private readonly CompactLcp _lcp;
    private readonly bool _tagged;

    /// <exception cref="ArgumentException">A tagged table over a reader of other than 4 bytes.</exception>
    internal ChildCells(TPos reader, CompactLcp lcp, bool tagged)
    {
        ArgumentNullException.ThrowIfNull(lcp);
        if (tagged && TPos.ElementBytes != 4) throw new ArgumentException("the tagged codec (CHILD codec 2) exists over 4-byte cells only (ADR-005 section 6.1)", nameof(tagged));
        _reader = reader;
        _lcp = lcp;
        _tagged = tagged;
    }

    /// <summary>Number of cells (<c>n + 1</c>).</summary>
    internal long Length => _reader.Count;

    /// <summary>True for CHILD codec 2.</summary>
    internal bool IsTagged => _tagged;

    /// <summary>The first l-index of the interval <c>[lb..rb]</c>: the cell at its left bound for a last child, at its right bound otherwise.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal long FirstLIndex(long lb, long rb, bool lastChild)
    {
        ulong cell = _reader.Get(lastChild ? lb : rb);
        return _tagged ? (long)(cell & TaggedValueMask) : (long)cell;
    }

    /// <summary>The l-index after <paramref name="k"/> in the interval that ends at <paramref name="rb"/>, or -1 if <paramref name="k"/> is the last one.</summary>
    internal long NextLIndex(long k, long rb)
    {
        if (_tagged)
        {
            ulong cell = _reader.Get(k);
            return (cell & TaggedKindMask) == TaggedNext ? (long)(cell & TaggedValueMask) : -1;
        }
        if (k >= rb) return -1;
        long lk = _lcp[k];
        if (_lcp[k + 1] == lk) return k + 1;
        long v = (long)_reader.Get(k);
        return _lcp[v] == lk ? v : -1;
    }
}
