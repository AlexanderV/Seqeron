using System.Runtime.CompilerServices;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa;

/// <summary>The sections that the navigation pack reads where they lie (work package S2.5): PSI and the persisted RMQ pyramid, each absent when the container has none.</summary>
internal sealed class NavSections(SaView? psi, SaView? rmq, int fanout, bool isChecked = true, long inPlaceFrom = EsaNavigationPack.InPlaceFrom)
{
    /// <summary>No section: an index built in memory.</summary>
    internal static NavSections None { get; } = new(null, null, SesaFormat.RmqFanout);

    /// <summary>The text length from which the 32-bit navigator of a container reads PSI in place instead of building the inverse suffix array (<see cref="EsaNavigationPack.InPlaceFrom"/>).</summary>
    internal long InPlaceFrom => inPlaceFrom;

    /// <summary>True in <see cref="EsaAccessMode.Checked"/> mode: what the pack reads from PSI is verified against the suffix array before it is used.</summary>
    internal bool IsChecked => isChecked;

    /// <summary>PSI: <c>PSI[i] = ISA[SA[i] + 1]</c> (0 for the empty suffix), <c>n + 1</c> entries.</summary>
    internal SaView? Psi => psi;

    /// <summary>The levels of the pyramid of the LCP array, level 0 first.</summary>
    internal SaView? Rmq => rmq;

    /// <summary>The fanout of the pyramid in <see cref="Rmq"/>.</summary>
    internal int Fanout => fanout;
}

/// <summary>
/// The navigation pack of an <see cref="EsaIndex"/> (ADR-001 R-SL, work packages 1b.4 and S2.5): what the suffix link of an internal node, <c>Parent</c>, <c>Lcp</c> and <c>Lca</c> need besides the suffix array and the LCP array. Kernel (c2) of
/// the Log of 1b.4a: for the node <c>v = [l..r]</c> of string depth <c>d</c> the link is the lcp-interval of value <c>t = d - 1</c> that contains <c>Psi(l)</c> and <c>Psi(r)</c> (fact F1), whose borders are <c>PrevLess(Psi(l), t)</c> and
/// <c>NextLess(Psi(r) + 1, t) - 1</c>; <c>Psi(i) = ISA[SA[i] + 1]</c>. The target's parent depth and last-child flag come out of the two LCP values just outside the target, <c>Pos = SA[Lb]</c> is one more read; there is no child-table walk, so
/// the kernel has no unbounded loop.
/// <para>Where the numbers come from (S2.5). <b>Psi</b> is the PSI section of a container, read where it lies (the builder writes it by default: no heap, no pass), else <c>ISA[SA[i] + 1]</c> from the inverse suffix array (built in memory:
/// eagerly for an index in memory, as in S1; on the first use for a container, and only for a text of at most <see cref="PackHeapLimit"/> symbols). <b>The pyramid</b> is the RMQ section of a container read where it lies, else built in memory
/// (<see cref="LcpPyramid"/>, 0.13 bytes per entry; <see cref="WideLcpPyramid"/> for the LCP arrays that the 32-bit form cannot hold). <c>Parent</c>, <c>SuffixLink</c>, <c>Lcp</c> and <c>Lca</c> of a container that has PSI and RMQ therefore
/// need no heap at all; <c>RankOf</c> and <c>ExtendLeft</c> need the inverse suffix array, so above the limit they are <see cref="NotSupportedException"/>, not an out-of-memory error.</para>
/// <para>Two kernels, one pack. The 32-bit kernel (<see cref="SuffixLink(in EsaNode)"/>, <see cref="Sa"/>, <see cref="Pyramid"/>) is the code of S1 for the navigator of the generic algorithms, over an index that the 32-bit navigator serves
/// (<see cref="EsaIndex.IsNarrow"/>); the long kernel (<see cref="SuffixLink(in EsaNode64)"/>, <see cref="Parent"/>, <see cref="Lca"/>, <see cref="Lcp"/>, <see cref="Isa(long, string)"/>) serves every index at any size and width, through
/// <see cref="PosSource"/> and <see cref="ILcpRanges"/>. Immutable after construction except the inverse suffix array, which is built once, thread-safe; safe for concurrent readers.</para>
/// </summary>
internal sealed unsafe class EsaNavigationPack : IDisposable
{
    /// <summary>
    /// The text length above which a container without an ISA section does not get a heap inverse suffix array (work package S2.5, ADR-005 section 3.3 and R5): 2^28 symbols, 1 GB of ISA. Above it <c>RankOf</c> and <c>ExtendLeft</c> are
    /// <see cref="NotSupportedException"/> that names the member, the limit and the way out; at or below it the array is built when it is first needed, as before. An index built in memory is not limited by it.
    /// </summary>
    internal const long PackHeapLimit = 1L << 28;

    /// <summary>
    /// The text length from which the 32-bit navigator of a container reads PSI in place instead of building the inverse suffix array (work package S2.5): 2^25 symbols, where the array would take 134 MB. Below it the pack is the pack of S1 (heap
    /// ISA, as the gate rows measured it): the walk of MEM, MUM and the anchors follows the query along the text, so its reads of ISA (indexed by text position) are consecutive, while PSI is indexed by rank and its reads are not; measured at 10^7 symbols
    /// the walk over PSI in place was 13 - 15 % slower than over the heap array (the log of S2.5). Above it the heap array is 4 bytes per symbol that a container does not need, and the walk reads PSI. The long kernel reads PSI at any size.
    /// </summary>
    internal const long InPlaceFrom = 1L << 25;

    private readonly EsaIndex? _index;
    private readonly long _n;
    private readonly bool _mapped;
    private readonly PosSource _sa64;
    private readonly PosSource? _psi;
    private readonly ILcpRanges _ranges;
    private readonly object _isaLock = new();
    private PinnedArray<int>? _isa;

    // the 32-bit kernel (an index that the 32-bit navigator serves)
    private readonly bool _narrow;
    private readonly int _n32;
    private readonly PinnedArray<int>? _sa;
    private readonly Buf<int> _saBuf;
    private readonly Buf<int> _isaBuf;
    private readonly Buf<int> _psiBuf;
    private readonly bool _viaPsi;
    private readonly bool _verifyPsi;
    private readonly bool _linkable;
    private readonly LcpPyramid? _pyramid;

    /// <summary>The pack of an index built in memory: the code of S1.</summary>
    private EsaNavigationPack(EsaIndex index, PinnedArray<int> sa, PinnedArray<int> isa, LcpPyramid pyramid)
    {
        _index = index;
        _n = index.TextLength;
        _n32 = checked((int)_n);
        _mapped = false;
        _narrow = true;
        _sa64 = PosSource.Of(index.Sa);
        _sa = sa;
        _saBuf = sa.Buf;
        _isa = isa;
        _isaBuf = isa.Buf;
        _linkable = true;
        _pyramid = pyramid;
        _ranges = pyramid;
    }

    /// <summary>The pack of a container (or of views of a test) that the 32-bit navigator serves.</summary>
    private EsaNavigationPack(EsaIndex index, Buf<int> sa, PinnedArray<int>? isa, Buf<int> psi, bool viaPsi, LcpPyramid pyramid)
    {
        _index = index;
        _n = index.TextLength;
        _n32 = checked((int)_n);
        _mapped = true;
        _narrow = true;
        _sa64 = PosSource.Of(index.Sa);
        _saBuf = sa;
        _isa = isa;
        _isaBuf = isa?.Buf ?? default;
        _psiBuf = psi;
        _viaPsi = viaPsi;
        _linkable = viaPsi || isa is not null;
        if (viaPsi)
        {
            NavSections sections = index.NavigationSections();
            _psi = PosSource.Of(sections.Psi!.Value);
            _verifyPsi = sections.IsChecked;
        }
        _pyramid = pyramid;
        _ranges = pyramid;
    }

    /// <summary>The pack of an index that only the long navigator serves, or that a test made over sources (<paramref name="index"/> null).</summary>
    private EsaNavigationPack(EsaIndex? index, long n, bool mapped, PosSource sa, PosSource? psi, ILcpRanges ranges, bool verifyPsi)
    {
        _verifyPsi = verifyPsi && psi is not null;
        _index = index;
        _n = n;
        _mapped = mapped;
        _narrow = false;
        _sa64 = sa;
        _psi = psi;
        _ranges = ranges;
    }

    /// <summary>Builds the pack of <paramref name="index"/> (O(n) for an index in memory: one scatter pass for ISA, one pass for the pyramid; O(1) plus what a container lacks for a container).</summary>
    /// <exception cref="ArgumentNullException"><paramref name="index"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The fanout is not a power of two of at least 2.</exception>
    /// <exception cref="InvalidDataException">A container whose PSI, RMQ or suffix array is damaged.</exception>
    internal static EsaNavigationPack Build(EsaIndex index, int fanout = LcpPyramid.DefaultFanout)
    {
        ArgumentNullException.ThrowIfNull(index);
        if (!index.IsMapped) return BuildMemory(index, fanout);
        NavSections sections = index.NavigationSections();
        int pyramidFanout = sections.Rmq is null ? fanout : sections.Fanout;
        if (!index.IsNarrow)
        {
            ILcpRanges ranges = LcpPyramid.Create(index.Lcp, pyramidFanout, sections.Rmq is { } rmq ? PosSource.Of(rmq) : null);
            return new EsaNavigationPack(index, index.TextLength, true, PosSource.Of(index.Sa), sections.Psi is { } psi ? PosSource.Of(psi) : null, ranges, sections.IsChecked);
        }
        return BuildNarrowMapped(index, sections, pyramidFanout);
    }

    /// <summary>The pack over sources that a test made (generated readers: no container, no <see cref="EsaIndex"/>): the long kernel only; the inverse suffix array does not exist.</summary>
    internal static EsaNavigationPack FromSources(long n, PosSource sa, PosSource? psi, ILcpRanges ranges, bool mapped = true, bool verifyPsi = true)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        ArgumentNullException.ThrowIfNull(sa);
        ArgumentNullException.ThrowIfNull(ranges);
        return new EsaNavigationPack(null, n, mapped, sa, psi, ranges, verifyPsi);
    }

    private static EsaNavigationPack BuildMemory(EsaIndex index, int fanout)
    {
        PinnedArray<int> sa = PinnedArray<int>.Pin(index.SuffixArray, "sa");
        PinnedArray<int>? isa = null;
        LcpPyramid? pyramid = null;
        try
        {
            pyramid = LcpPyramid.Build(index.Lcp, fanout);
            isa = PinnedArray<int>.Allocate(index.SuffixArray.Length, zero: false, "isa");
            Buf<int> saBuf = sa.Buf;
            Buf<int> isaBuf = isa.Buf;
            for (int rank = 0, count = index.SuffixArray.Length; rank < count; rank++) isaBuf[saBuf[rank]] = rank;
            return new EsaNavigationPack(index, sa, isa, pyramid);
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
    /// The pack of a container that the 32-bit navigator serves. The suffix array is read where it lies. With PSI the endpoints of a link are two reads of PSI and nothing is built; without it the inverse suffix array is built (the entries are
    /// untrusted: the scatter checks every one, a value outside <c>0..n</c> or one that occurs twice is an <see cref="InvalidDataException"/>, and after n + 1 distinct values every slot holds a rank, which makes the later reads of ISA safe), unless the
    /// text is above <see cref="PackHeapLimit"/> symbols, when the suffix link is <see cref="NotSupportedException"/> by name. The pyramid is the RMQ section where it lies, else it is built.
    /// </summary>
    private static EsaNavigationPack BuildNarrowMapped(EsaIndex index, NavSections sections, int fanout)
    {
        SaView sa = index.Sa;
        PinnedArray<int>? isa = null;
        LcpPyramid? pyramid = null;
        try
        {
            pyramid = sections.Rmq is { } rmq
                ? LcpPyramid.FromLevels(index.Lcp, fanout, rmq.Pointer, rmq.Length)
                : LcpPyramid.Build(index.Lcp, fanout);
            int count = checked((int)sa.Length); // a narrow container has at most int.MaxValue ranks
            Buf<int> psiBuf = default;
            bool viaPsi = sections.Psi is not null && (index.TextLength >= sections.InPlaceFrom || index.TextLength > PackHeapLimit);
            if (viaPsi) psiBuf = new Buf<int>(sections.Psi!.Value.Pointer, count, "psi");
            else if (index.TextLength <= PackHeapLimit) isa = BuildIsaOf(PosSource.Of(sa), count);
            return new EsaNavigationPack(index, new Buf<int>(sa.Pointer, count, "sa"), isa, psiBuf, viaPsi, pyramid);
        }
        catch
        {
            isa?.Dispose();
            pyramid?.Dispose();
            throw;
        }
    }

    /// <summary>The inverse of the suffix array <paramref name="sa"/> (<paramref name="count"/> = n + 1 entries), checked: the entries of a container are untrusted.</summary>
    private static PinnedArray<int> BuildIsaOf(PosSource sa, int count)
    {
        PinnedArray<int> isa = PinnedArray<int>.Allocate(count, zero: false, "isa");
        try
        {
            Buf<int> isaBuf = isa.Buf;
            isaBuf.Fill(-1, count);
            for (int rank = 0; rank < count; rank++)
            {
                long position = sa[rank];
                if ((ulong)position >= (ulong)count || isaBuf[checked((int)position)] != -1)
                    throw new InvalidDataException($"the suffix array holds {position} at rank {rank}: it is not a permutation of 0..{count - 1}, the file is damaged");
                isaBuf[checked((int)position)] = rank;
            }
            return isa;
        }
        catch
        {
            isa.Dispose();
            throw;
        }
    }

    /// <summary>The index the pack belongs to (null for a pack that a test made over sources).</summary>
    internal EsaIndex Index => _index ?? throw new NotSupportedException("this pack was made over sources and belongs to no index");

    /// <summary>The root node <c>(0, n, 0, 0, n, true)</c>.</summary>
    internal EsaNode Root => new(0, _n32, 0, 0, _n32, true);

    /// <summary>The root at 64 bits.</summary>
    internal EsaNode64 Root64 => new(0, _n, 0, 0, _n, true);

    /// <summary>Number of symbols n.</summary>
    internal long Length => _n;

    /// <summary>The LCP helper as the 64-bit contract has it (the 32-bit form or the wide one).</summary>
    internal ILcpRanges Ranges => _ranges;

    /// <summary>The LCP helper of the 32-bit navigator.</summary>
    /// <exception cref="NotSupportedException">The index is not served by the 32-bit navigator.</exception>
    internal LcpPyramid Pyramid => _pyramid ?? throw new NotSupportedException("the 32-bit pyramid exists for the indexes that the 32-bit navigator serves; use Ranges");

    /// <summary>The suffix array at 64-bit ranks.</summary>
    internal PosSource Sa64 => _sa64;

    /// <summary>The pinned suffix array (<c>n + 1</c> entries) that the 32-bit navigator reads.</summary>
    internal Buf<int> Sa => _narrow ? _saBuf : throw new NotSupportedException("the 32-bit suffix array view exists for the indexes that the 32-bit navigator serves; use Sa64");

    /// <summary>
    /// True when the 32-bit suffix link can run: the pack has PSI where it lies or an inverse suffix array. False only for a container without PSI of more than <see cref="PackHeapLimit"/> symbols (the 32-bit navigator, which cannot link,
    /// is refused at its creation by <see cref="EnsureLinkable"/>; the long kernel says the same thing by name at the first link).
    /// </summary>
    internal bool IsLinkable => _linkable || !_narrow;

    /// <summary>Refuses a navigator that could not link.</summary>
    /// <exception cref="NotSupportedException">The pack has neither PSI nor an inverse suffix array (see <see cref="IsLinkable"/>).</exception>
    internal void EnsureLinkable()
    {
        if (_narrow && !_linkable) throw IsaMissing("the suffix link (MEM, MUM, anchors and the longest common substrings walk by suffix links)");
    }

    /// <summary>True when the suffix link reads the PSI section where it lies (no inverse suffix array is needed for it).</summary>
    internal bool UsesPsi => _psi is not null;

    /// <summary>True once the inverse suffix array exists (always for an index built in memory; for a container after the first <c>RankOf</c> or <c>ExtendLeft</c>, or a suffix link without PSI).</summary>
    internal bool IsaBuilt => Volatile.Read(ref _isa) is not null;

    /// <summary><c>ISA[position]</c> (the 32-bit navigator: the rank of the suffix that starts at <paramref name="position"/>, <c>0..n</c>; <c>ISA[n] = 0</c>).</summary>
    internal int Isa(int position)
    {
        if ((uint)position > (uint)_n32)
        {
            if (_mapped) throw new InvalidDataException($"the inverse suffix array is read at position {position}, but the text has {_n} symbols: the file is damaged");
            throw new ArgumentOutOfRangeException(nameof(position));
        }
        return _isaBuf[position];
    }

    /// <summary><c>ISA[position]</c> at 64 bits: the inverse suffix array is built on the first call (below <see cref="PackHeapLimit"/> for a container).</summary>
    /// <exception cref="NotSupportedException">A container of more than <see cref="PackHeapLimit"/> symbols: the member, the limit and the way out are named.</exception>
    internal long Isa(long position, string member)
    {
        if ((ulong)position > (ulong)_n)
        {
            if (_mapped) throw new InvalidDataException($"the inverse suffix array is read at position {position}, but the text has {_n} symbols: the file is damaged");
            throw new ArgumentOutOfRangeException(nameof(position));
        }
        return GetIsa(member).Buf[checked((int)position)];
    }

    /// <summary>The inverse suffix array as a managed array (tests and diagnostics); built when it does not exist.</summary>
    internal int[] IsaArray => GetIsa("EsaNavigationPack.IsaArray").Array;

    private PinnedArray<int> GetIsa(string member)
    {
        PinnedArray<int>? isa = Volatile.Read(ref _isa);
        if (isa is not null) return isa;
        if (_mapped && _n > PackHeapLimit)
            throw new NotSupportedException($"{member} needs the inverse suffix array (ISA) of the text: a container of {_n} symbols has no ISA section (this version of the format writes none) and the heap holds ISA only for a text of at most PackHeapLimit = {PackHeapLimit} symbols (4 bytes per symbol, 1 GB). The suffix link, Parent, Lcp, Lca and SuffixAt need only the PSI and RMQ sections that the builder writes, and work at any size; for RankOf and ExtendLeft use a text of at most {PackHeapLimit} symbols");
        if (_index is null) throw new NotSupportedException($"{member}: this pack was made over sources and has no inverse suffix array");
        lock (_isaLock)
        {
            isa = _isa;
            if (isa is null)
            {
                isa = BuildIsaOf(_sa64, checked((int)(_n + 1)));
                Volatile.Write(ref _isa, isa);
            }
            return isa;
        }
    }

    /// <summary>Bytes that the pack holds on the heap: the inverse suffix array if it exists, plus the pyramid levels that were built (0 for the levels of the RMQ section).</summary>
    internal long PayloadBytes => (4L * (Volatile.Read(ref _isa)?.Length ?? 0)) + _ranges.LevelBytes;

    // ------------------------------------------------------------------------------------------------ the 32-bit kernel (the code of S1, with PSI as a second source of the endpoints)

    /// <summary>
    /// The suffix link of <paramref name="node"/>: for an internal node <c>aY</c> the node <c>Y</c> (the root for a node of depth 1), for a leaf the leaf of the next suffix (the root for the terminator-only leaf), for the root the root. The
    /// result carries its parent depth, last-child flag and <c>Pos</c>.
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
        int a, b;
        if (_viaPsi)
        {
            if (_mapped && ((uint)node.Rb > (uint)_n32 || (uint)node.Lb > (uint)_n32)) throw DamagedNode();
            a = _psiBuf[node.Lb];
            b = _psiBuf[node.Rb];
            // PSI is untrusted: in Checked mode an endpoint must be the rank of the next suffix (SA[Psi(i)] = SA[i] + 1), or the link would be a node of the tree that is not the link; a rank outside 0..n is damage in both modes
            if (_mapped && ((uint)a > (uint)_n32 || (uint)b > (uint)_n32)) throw DamagedNode();
            if (_verifyPsi && (sa[a] != node.Pos + 1 || sa[b] != sa[node.Rb] + 1)) throw DamagedNode();
        }
        else
        {
            Buf<int> isa = _isaBuf;
            if (_mapped && ((uint)node.Rb > (uint)_n32 || (uint)(node.Pos + 1) > (uint)_n32 || (uint)(sa[node.Rb] + 1) > (uint)_n32)) throw DamagedNode();
            a = isa[node.Pos + 1];
            b = isa[sa[node.Rb] + 1];
        }
        if (_mapped && (b < a || (uint)a > (uint)_n32 || (uint)b > (uint)_n32 || (uint)node.Lb > (uint)_n32)) throw DamagedNode();
        _pyramid!.Borders(a, b, t, out int left, out int right, out int leftValue, out int rightValue);
        int parentDepth = Math.Max(leftValue, rightValue);
        if (_mapped && (parentDepth >= t || (uint)left > (uint)_n32)) throw DamagedNode(); // the target of a link has a longer edge than a point: otherwise the walks of the generic algorithms would not advance
        int target = sa[left];
        if (_mapped && (uint)target > (uint)_n32) throw DamagedNode();
        return new EsaNode(left, right - 1, t, parentDepth, target, leftValue > rightValue);
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private EsaNode LeafLink(in EsaNode leaf)
    {
        int position = leaf.Pos;
        if (position >= _n32) return Root; // the terminator-only leaf of the root
        if (_mapped && position < 0) throw DamagedNode();
        int rank;
        if (_viaPsi)
        {
            if (_mapped && (uint)leaf.Lb > (uint)_n32) throw DamagedNode();
            rank = _psiBuf[leaf.Lb];
            if (_mapped && (uint)rank > (uint)_n32) throw DamagedNode();
            if (_verifyPsi && _saBuf[rank] != position + 1) throw DamagedNode();
        }
        else
        {
            rank = _isaBuf[position + 1];
        }
        int left = _pyramid!.ValueAt(rank);
        int right = rank < _n32 ? _pyramid.ValueAt(rank + 1) : -1;
        if (_mapped && Math.Max(left, right) >= leaf.Depth - 1) throw DamagedNode();
        return new EsaNode(rank, rank, leaf.Depth - 1, Math.Max(left, right), position + 1, left > right);
    }

    // ------------------------------------------------------------------------------------------------ the long kernel (any size, any width)

    /// <summary><c>Psi(rank) = ISA[SA[rank] + 1]</c>, 0 for the empty suffix: from the PSI section where it lies, else from the inverse suffix array.</summary>
    internal long Psi(long rank, string member)
    {
        if ((ulong)rank > (ulong)_n) throw _mapped ? DamagedNode() : new ArgumentOutOfRangeException(nameof(rank));
        if (_psi is not null) return _psi[rank];
        long position = _sa64[rank];
        if (position >= _n) return 0;
        if (position < 0 || position + 1 > _n) throw DamagedNode();
        return Isa(position + 1, member);
    }

    /// <summary>
    /// The suffix link of <paramref name="node"/> at 64 bits (see <see cref="SuffixLink(in EsaNode)"/>): the endpoints are PSI of the two bounds (PSI where it lies, else ISA of <c>SA[i] + 1</c>), the borders come from the pyramid at
    /// 64-bit ranks and values.
    /// </summary>
    /// <exception cref="ArgumentException"><paramref name="node"/> is the null node.</exception>
    /// <exception cref="NotSupportedException">A container without PSI of more than <see cref="PackHeapLimit"/> symbols: the link needs the inverse suffix array.</exception>
    /// <exception cref="InvalidDataException">The container is damaged where the link reads it.</exception>
    internal EsaNode64 SuffixLink(in EsaNode64 node)
    {
        if (node.Lb < 0) throw new ArgumentException("the null node has no suffix link", nameof(node));
        if (node.Depth == 0) return Root64;
        if (node.Lb == node.Rb) return LeafLink(node);
        long t = node.Depth - 1;
        if (t == 0) return Root64;
        long a = Psi(node.Lb, "the suffix link");
        long b = Psi(node.Rb, "the suffix link");
        if (_mapped && (a > _n || b < a || b > _n)) throw DamagedNode();
        if (_verifyPsi && (_sa64[a] != node.Pos + 1 || _sa64[b] != _sa64[node.Rb] + 1)) throw DamagedNode();
        _ranges.Borders(a, b, t, out long left, out long right, out long leftValue, out long rightValue);
        long parentDepth = Math.Max(leftValue, rightValue);
        if (_mapped && (parentDepth >= t || left < 0 || left > _n)) throw DamagedNode();
        long target = _sa64[left];
        if (_mapped && (ulong)target > (ulong)_n) throw DamagedNode();
        return new EsaNode64(left, right - 1, t, parentDepth, target, leftValue > rightValue);
    }

    private EsaNode64 LeafLink(in EsaNode64 leaf)
    {
        long position = leaf.Pos;
        if (position >= _n) return Root64;
        if (_mapped && position < 0) throw DamagedNode();
        long rank = Psi(leaf.Lb, "the suffix link of a leaf");
        if (_mapped && rank > _n) throw DamagedNode();
        if (_verifyPsi && _sa64[rank] != position + 1) throw DamagedNode();
        long left = _ranges.ValueAt(rank);
        long right = rank < _n ? _ranges.ValueAt(rank + 1) : -1;
        if (_mapped && Math.Max(left, right) >= leaf.Depth - 1) throw DamagedNode();
        return new EsaNode64(rank, rank, leaf.Depth - 1, Math.Max(left, right), position + 1, left > right);
    }

    /// <summary>
    /// The nearest explicit ancestor of an explicit node (the root for a node whose parent depth is 0): the lcp-interval of value <c>ParentDepth</c> around it, whose borders are the previous and the next LCP value below that depth.
    /// </summary>
    /// <exception cref="ArgumentException"><paramref name="node"/> is the root (it has no parent) or the null node.</exception>
    internal EsaNode64 Parent(in EsaNode64 node)
    {
        if (node.Lb < 0 || node.Depth == 0) throw new ArgumentException("the root and the null node have no parent", nameof(node));
        long parentDepth = node.ParentDepth;
        if (parentDepth == 0) return Root64;
        _ranges.Borders(node.Lb, node.Rb, parentDepth, out long left, out long right, out long leftValue, out long rightValue);
        return new EsaNode64(left, right - 1, parentDepth, Math.Max(leftValue, rightValue), left >= 0 && left <= _n ? _sa64[left] : -1, leftValue > rightValue);
    }

    /// <summary><c>min LCP[min + 1 .. max]</c> of two ranks (the Lcp of the two suffixes); equal ranks give the length of the suffix.</summary>
    internal long Lcp(long rankA, long rankB)
    {
        if (rankA == rankB) return _n - _sa64[rankA];
        return _ranges.Min(Math.Min(rankA, rankB) + 1, Math.Max(rankA, rankB));
    }

    /// <summary>
    /// The lowest common ancestor of the leaves of two ranks: the lcp-interval of value <c>m = min LCP[lo + 1 .. hi]</c> around them, whose borders are the previous and the next LCP value below <c>m</c> (the root when <c>m = 0</c>); the
    /// leaf itself for equal ranks, the root for the empty text. The result is <c>(Lb, Rb, Depth)</c>.
    /// </summary>
    internal (long Lb, long Rb, long Depth) Lca(long rankA, long rankB)
    {
        if (_n == 0) return (0, 0, 0);
        if (rankA == rankB)
        {
            long depth = _n - _sa64[rankA] + 1;
            return (rankA, rankA, depth);
        }
        long lo = Math.Min(rankA, rankB);
        long hi = Math.Max(rankA, rankB);
        long value = _ranges.Min(lo + 1, hi);
        if (value == 0) return (0, _n, 0);
        return (_ranges.PrevLess(lo, value), _ranges.NextLess(hi + 1, value) - 1, value);
    }

    /// <summary>Releases the pins; idempotent. The pack must not be used afterwards.</summary>
    public void Dispose()
    {
        _sa?.Dispose();
        _isa?.Dispose();
        _ranges.Dispose();
    }

    private static InvalidDataException DamagedNode() => new("a node of the container's tree has no suffix link in the inverse suffix array or in PSI: the file is damaged");

    private NotSupportedException IsaMissing(string use) =>
        new($"{use} needs PSI or the inverse suffix array: this container has no PSI section and has {_n} symbols, more than the PackHeapLimit = {PackHeapLimit} of the heap array (4 bytes per symbol). Build the container with PSI (the default of the builder)");
}
