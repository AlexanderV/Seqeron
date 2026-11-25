using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

/// <summary>
/// The navigation pack of an <see cref="EsaIndex"/> (ADR-001 R-SL, work package 1b.4): the inverse suffix array (4 bytes per rank) and the <see cref="LcpPyramid"/> (0.13 bytes per rank at
/// fanout 32), which together give the suffix link of an internal node. Kernel (c2) of the Log of 1b.4a: for the node <c>v = [l..r]</c> of string depth <c>d</c> the link is the lcp-interval of
/// value <c>t = d - 1</c> that contains <c>Psi(l)</c> and <c>Psi(r)</c> (fact F1), whose borders are <c>PrevLess(Psi(l), t)</c> and <c>NextLess(Psi(r) + 1, t) - 1</c>; <c>Psi(i) = ISA[SA[i] + 1]</c>,
/// so the endpoints are read through ISA (no Psi array: measured 6-24 % faster end to end). The target's parent depth and last-child flag come out of the two LCP values just outside the target,
/// <c>Pos = SA[Lb]</c> is one more read; there is no child-table walk, so the kernel has no unbounded loop. Built by <see cref="EsaIndex.Navigation"/> (lazily, once, thread-safe, heap-resident); it
/// also holds the pinned view of the suffix array that the navigator reads through pointers. Immutable; safe for concurrent readers.
/// </summary>
internal sealed unsafe class EsaNavigationPack : IDisposable
{
    /// <summary>
    /// The text length above which a container without an ISA section does not get a heap inverse suffix array (work package S2.5, ADR-005 section 3.3 and R5): 2^28 symbols, 1 GB of ISA. Above it <c>RankOf</c> and <c>ExtendLeft</c> are
    /// <see cref="NotSupportedException"/> that names the member, the limit and the way out; at or below it the array is built when it is first needed, as before.
    /// </summary>
    internal const long PackHeapLimit = 1L << 28;

    private readonly EsaIndex _index;
    private readonly int _n;
    private readonly PinnedArray<int>? _sa;
    private readonly PinnedArray<int> _isa;
    private readonly bool _mapped;
    private readonly Buf<int> _saBuf;
    private readonly Buf<int> _isaBuf;
    private readonly LcpPyramid _pyramid;

    private EsaNavigationPack(EsaIndex index, PinnedArray<int>? sa, Buf<int> saBuf, PinnedArray<int> isa, LcpPyramid pyramid)
    {
        _index = index;
        _n = index.Length;
        _sa = sa;
        _mapped = index.IsMapped;
        _isa = isa;
        _saBuf = saBuf;
        _isaBuf = isa.Buf;
        _pyramid = pyramid;
    }

    /// <summary>Builds the pack of <paramref name="index"/> (O(n): one scatter pass for ISA, one pass for the pyramid).</summary>
    /// <exception cref="ArgumentNullException"><paramref name="index"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The fanout is not a power of two of at least 2.</exception>
    internal static EsaNavigationPack Build(EsaIndex index, int fanout = LcpPyramid.DefaultFanout)
    {
        ArgumentNullException.ThrowIfNull(index);
        if (index.IsMapped) return BuildMapped(index, fanout);
        var sa = PinnedArray<int>.Pin(index.SuffixArray, "sa");
        PinnedArray<int>? isa = null;
        LcpPyramid? pyramid = null;
        try
        {
            pyramid = LcpPyramid.Build(index.Lcp, fanout);
            isa = PinnedArray<int>.Allocate(index.SuffixArray.Length, zero: false, "isa");
            Buf<int> saBuf = sa.Buf;
            Buf<int> isaBuf = isa.Buf;
            for (int rank = 0, count = index.SuffixArray.Length; rank < count; rank++) isaBuf[saBuf[rank]] = rank;
            return new EsaNavigationPack(index, sa, saBuf, isa, pyramid);
        }
        catch
        {
            sa.Dispose();
            isa?.Dispose();
            pyramid?.Dispose();
            throw;
        }
    }

    /// <summary>
    /// The pack of an index that reads a container: the suffix array is read where it lies (nothing is pinned or copied), the inverse is built in memory like the in-memory index's. The entries are untrusted, so the scatter
    /// checks every one: a value outside <c>0..n</c>, or one that occurs twice, is an <see cref="InvalidDataException"/> (the slots start at -1, so after n + 1 distinct in-range values every slot holds a rank, which is what makes
    /// the later reads of ISA safe).
    /// </summary>
    private static EsaNavigationPack BuildMapped(EsaIndex index, int fanout)
    {
        SaView sa = index.Sa;
        PinnedArray<int>? isa = null;
        LcpPyramid? pyramid = null;
        try
        {
            pyramid = LcpPyramid.Build(index.Lcp, fanout);
            int count = IntBridge.Length(sa.Length, "the suffix array of the container");
            isa = PinnedArray<int>.Allocate(count, zero: false, "isa");
            Buf<int> isaBuf = isa.Buf;
            isaBuf.Fill(-1, count);
            for (int rank = 0; rank < count; rank++)
            {
                int position = sa[rank];
                if ((uint)position >= (uint)count || isaBuf[position] != -1)
                    throw new InvalidDataException($"the suffix array holds {position} at rank {rank}: it is not a permutation of 0..{count - 1}, the file is damaged");
                isaBuf[position] = rank;
            }
            return new EsaNavigationPack(index, null, new Buf<int>(sa.Pointer, count, "sa"), isa, pyramid);
        }
        catch
        {
            isa?.Dispose();
            pyramid?.Dispose();
            throw;
        }
    }

    /// <summary>The index the pack belongs to.</summary>
    internal EsaIndex Index => _index;

    /// <summary>The root node <c>(0, n, 0, 0, n, true)</c>.</summary>
    internal EsaNode Root => new(0, _n, 0, 0, _n, true);

    /// <summary>The LCP helper.</summary>
    internal LcpPyramid Pyramid => _pyramid;

    /// <summary>The pinned suffix array (<c>n + 1</c> entries).</summary>
    internal Buf<int> Sa => _saBuf;

    /// <summary><c>ISA[position]</c>, the rank of the suffix that starts at <paramref name="position"/> (<c>0..n</c>; <c>ISA[n] = 0</c>).</summary>
    internal int Isa(int position)
    {
        if ((uint)position > (uint)_n)
        {
            if (_mapped) throw new InvalidDataException($"the inverse suffix array is read at position {position}, but the text has {_n} symbols: the file is damaged");
            throw new ArgumentOutOfRangeException(nameof(position));
        }
        return _isaBuf[position];
    }

    /// <summary>The inverse suffix array as a managed array (tests and diagnostics).</summary>
    internal int[] IsaArray => _isa.Array;

    /// <summary>Bytes of ISA plus the pyramid levels (what the pack adds to the index).</summary>
    internal long PayloadBytes => (4L * _isa.Length) + _pyramid.LevelBytes;

    /// <summary>
    /// The suffix link of <paramref name="node"/>: for an internal node <c>aY</c> the node <c>Y</c> (the root for a node of depth 1), for a leaf the leaf of the next suffix (the root for the
    /// terminator-only leaf), for the root the root. The result carries its parent depth, last-child flag and <c>Pos</c>.
    /// </summary>
    /// <exception cref="ArgumentException"><paramref name="node"/> is the null node.</exception>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal EsaNode SuffixLink(in EsaNode node)
    {
        if (node.Lb < 0) throw new ArgumentException("the null node has no suffix link", nameof(node));
        if (node.Depth == 0) return Root;                        // the root
        if (node.Lb == node.Rb) return LeafLink(node);           // a leaf: the leaf of the next suffix
        int t = node.Depth - 1;
        if (t == 0) return Root;                                  // an internal node of depth 1

        Buf<int> sa = _saBuf;
        Buf<int> isa = _isaBuf;
        if (_mapped && ((uint)node.Rb > (uint)_n || (uint)(node.Pos + 1) > (uint)_n || (uint)(sa[node.Rb] + 1) > (uint)_n)) throw DamagedNode();
        int a = isa[node.Pos + 1];
        int b = isa[sa[node.Rb] + 1];
        if (_mapped && (b < a || (uint)node.Lb > (uint)_n)) throw DamagedNode();
        _pyramid.Borders(a, b, t, out int left, out int right, out int leftValue, out int rightValue);
        int parentDepth = Math.Max(leftValue, rightValue);
        if (_mapped && (parentDepth >= t || (uint)left > (uint)_n)) throw DamagedNode(); // the target of a link has a longer edge than a point: otherwise the walks of the generic algorithms would not advance
        return new EsaNode(left, right - 1, t, parentDepth, sa[left], leftValue > rightValue);
    }

    /// <summary>Releases the pins; idempotent. The pack must not be used afterwards.</summary>
    public void Dispose()
    {
        _sa?.Dispose();
        _isa.Dispose();
        _pyramid.Dispose();
    }

    private static InvalidDataException DamagedNode() => new("a node of the container's tree has no suffix link in the inverse suffix array: the file is damaged");

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private EsaNode LeafLink(in EsaNode leaf)
    {
        int position = leaf.Pos;
        if (position >= _n) return Root; // the terminator-only leaf of the root
        if (_mapped && position < 0) throw DamagedNode();
        int rank = _isaBuf[position + 1];
        int left = _pyramid.ValueAt(rank);
        int right = rank < _n ? _pyramid.ValueAt(rank + 1) : -1;
        if (_mapped && Math.Max(left, right) >= leaf.Depth - 1) throw DamagedNode();
        return new EsaNode(rank, rank, leaf.Depth - 1, Math.Max(left, right), position + 1, left > right);
    }
}
