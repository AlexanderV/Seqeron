using System.Runtime.CompilerServices;

namespace SuffixTree.Esa.Bench.Links;

/// <summary>Counters of one <see cref="MemDriver.Run"/>: query starts, suffix links followed (kernel calls), successful extension steps (one query symbol consumed each), child lookups made to extend from a node, child lookups made by the rescan after a link, matches emitted.</summary>
internal sealed record MemStats(long Starts, long Links, long ExtendSteps, long ChildLookups, long RescanLookups, long Emitted);

/// <summary>
/// The throwaway matching-statistics driver of prototype P2 (work package 1b.4a, not the navigator of 1b.4): MEM enumeration on the ESA with the structure of
/// <c>SuffixTreeAlgorithms.FindMaximalExactMatches</c> (extend the match, follow the suffix link of the last node, rescan down to the new locus) where a locus is a rank
/// interval and the suffix link comes from a <see cref="LinkKernel"/>. A match is reported by scanning the LCP array outward from the rank interval of the longest match
/// while the running minimum stays at least <c>minLength</c>: those are exactly the suffixes whose common prefix with the query reaches <c>minLength</c>.
/// </summary>
internal sealed class MemDriver
{
    private struct Node
    {
        internal int Lb;
        internal int Rb;
        internal int Depth;
        internal int Pos;
        internal bool Last;
    }

    private struct Locus
    {
        /// <summary>The deepest node on the path with depth at most <see cref="Len"/> (the root for the empty match).</summary>
        internal Node V;

        /// <summary>The child of <see cref="V"/> whose edge holds the match, valid when <see cref="HasEdge"/> (then <c>V.Depth &lt; Len &lt; C.Depth</c>).</summary>
        internal Node C;
        internal bool HasEdge;
        internal int Len;
    }

    private readonly LinkContext _links;
    private readonly LinkFunc<NoCount> _link;
    private readonly int _loci;
    private readonly int _n;
    private readonly int[] _sa;
    private readonly byte[] _text;
    private readonly CompactLcp _lcp;
    private readonly TaggedChild _child;
    private readonly int[] _symbolOf;
    private long _starts, _followed, _extend, _lookups, _rescan, _emitted;
    private long _linkTicks;

    /// <param name="index">The ESA (tagged child table, at most 256 distinct symbols).</param>
    /// <param name="links">ISA / Psi / helper for the kernel.</param>
    /// <param name="kernel">The link kernel.</param>
    /// <param name="source">The endpoint source.</param>
    /// <param name="loci">1: one chain (matching statistics of the ESA); 2: the generic algorithm's two loci, the second chain only costs its suffix links (cost emulation).</param>
    internal MemDriver(EsaIndex index, LinkContext links, LinkKernel kernel, EndpointSource source, int loci)
    {
        ArgumentNullException.ThrowIfNull(index);
        ArgumentNullException.ThrowIfNull(links);
        if (loci is not (1 or 2)) throw new ArgumentOutOfRangeException(nameof(loci), loci, "1 or 2 loci");
        if (index.Prepared.Width != SymbolWidth.Byte) throw new NotSupportedException("the driver reads byte symbols (at most 256 distinct code units)");
        _links = links;
        _link = LinkKernels.Get<NoCount>(kernel, source);
        _loci = loci;
        _n = index.Length;
        _sa = index.SuffixArray;
        _text = index.Prepared.Bytes.ToArray();
        _lcp = index.Lcp;
        _child = links.Child;
        char[] alphabet = index.Prepared.Alphabet;
        _symbolOf = alphabet.Length == 0 ? [] : new int[alphabet[^1] + 1];
        Array.Fill(_symbolOf, -2);
        for (int s = 0; s < alphabet.Length; s++) _symbolOf[alphabet[s]] = s;
    }

    /// <summary>When true every kernel call is timed (two timestamps per link); the sum is <see cref="LinkSeconds"/>. Off by default.</summary>
    internal bool TimeLinks { get; set; }

    /// <summary>Seconds spent inside the link kernel during the last <see cref="Run"/> (only with <see cref="TimeLinks"/>).</summary>
    internal double LinkSeconds => _linkTicks / (double)System.Diagnostics.Stopwatch.Frequency;

    /// <summary>Counters of the last <see cref="Run"/>.</summary>
    internal MemStats Stats => new(_starts, _followed, _extend, _lookups, _rescan, _emitted);

    /// <summary>All maximal exact matches of at least <paramref name="minLength"/> symbols, ascending by query position, then text position.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> Run(string query, int minLength)
    {
        ArgumentNullException.ThrowIfNull(query);
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 1);
        _starts = _followed = _extend = _lookups = _rescan = _emitted = 0;
        _linkTicks = 0;
        var results = new List<(int PositionInText, int PositionInQuery, int Length)>();
        int q = query.Length;
        if (q < minLength || _n < minLength) return results;

        var qsym = new int[q];
        for (int i = 0; i < q; i++)
        {
            char c = query[i];
            qsym[i] = c < _symbolOf.Length ? _symbolOf[c] : -2;
        }

        Locus m = Root();
        Locus p = Root();
        var hits = new List<(int Position, int Length)>();
        for (int s = 0; s < q; s++)
        {
            _starts++;
            while (s + m.Len < q && TryExtend(ref m, qsym[s + m.Len])) _extend++;

            if (_loci == 2)
            {
                // the generic algorithm's second locus: the string of length min(ms, minLength)
                int target = Math.Min(m.Len, minLength);
                while (p.Len < target)
                {
                    if (!TryExtend(ref p, qsym[s + p.Len])) throw new InvalidOperationException("the prefix of a match did not extend");
                    _extend++;
                }
            }

            if (m.Len >= minLength) Emit(m, s, qsym, minLength, hits, results);

            if (m.Len > 0) FollowLink(ref m, s, qsym);
            if (_loci == 2 && p.Len > 0) FollowLink(ref p, s, qsym);
        }
        return results;
    }

    private Locus Root() => new() { V = new Node { Lb = 0, Rb = _n, Depth = 0, Pos = _sa[0], Last = true } };

    /// <summary>The child of <paramref name="v"/> whose edge starts with the dense symbol <paramref name="sym"/> (children ascend by key, the terminator first), as <c>Descend</c> of the search finds it.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private bool TryGetChild(in Node v, int sym, out Node child)
    {
        int a = v.Lb;
        int boundary = _child.FirstLIndex(v.Lb, v.Rb, v.Last);
        int b;
        while (true)
        {
            b = boundary >= 0 ? boundary - 1 : v.Rb;
            int p = _sa[a] + v.Depth;
            int key = p < _n ? _text[p] : -1;
            if (key == sym) break;
            if (key > sym || boundary < 0)
            {
                child = default;
                return false;
            }
            a = boundary;
            boundary = _child.NextLIndex(boundary, v.Rb);
        }

        bool last = boundary < 0;
        int pos = _sa[a];
        child = new Node
        {
            Lb = a,
            Rb = b,
            Pos = pos,
            Last = last,
            Depth = a == b ? _n - pos + 1 : _lcp[_child.FirstLIndex(a, b, last)],
        };
        return true;
    }

    /// <summary>Extends the match by one symbol (<paramref name="sym"/> &lt; 0: a symbol the text does not contain).</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private bool TryExtend(ref Locus l, int sym)
    {
        if (sym < 0) return false;
        if (l.HasEdge)
        {
            int p = l.C.Pos + l.Len;
            int key = p < _n ? _text[p] : -1;
            if (key != sym) return false;
            l.Len++;
            if (l.Len == l.C.Depth)
            {
                l.V = l.C;
                l.HasEdge = false;
            }
            return true;
        }

        if (!TryGetChild(l.V, sym, out Node c)) return false;
        _lookups++;
        l.C = c;
        l.HasEdge = true;
        l.Len++;
        if (l.Len == c.Depth)
        {
            l.V = c;
            l.HasEdge = false;
        }
        return true;
    }

    /// <summary>
    /// The locus of the match without its first symbol: the suffix link of the node above the locus, then the rescan down the path of the remaining symbols (skip/count: only
    /// the first symbol of every edge is looked at). <paramref name="start"/> is the query position where the old match began.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private void FollowLink(ref Locus l, int start, int[] qsym)
    {
        if (l.V.Depth > 0)
        {
            var node = new LinkNode(l.V.Lb, l.V.Rb, l.V.Depth, l.V.Pos);
            var none = new NoCount();
            LinkTarget t;
            if (TimeLinks)
            {
                long t0 = System.Diagnostics.Stopwatch.GetTimestamp();
                _link(_links, node, ref none, out t);
                _linkTicks += System.Diagnostics.Stopwatch.GetTimestamp() - t0;
            }
            else
            {
                _link(_links, node, ref none, out t);
            }
            _followed++;
            l.V = new Node { Lb = t.Lb, Rb = t.Rb, Depth = l.V.Depth - 1, Last = t.IsLast, Pos = _sa[t.Lb] };
        }
        l.Len--;
        l.HasEdge = false;

        int remaining = l.Len - l.V.Depth;
        int pos = start + 1 + l.V.Depth;
        while (remaining > 0)
        {
            if (!TryGetChild(l.V, qsym[pos], out Node c)) throw new InvalidOperationException("the rescan left the tree: a suffix of a substring must be a substring");
            _rescan++;
            int edge = c.Depth - l.V.Depth;
            if (edge <= remaining)
            {
                l.V = c;
                pos += edge;
                remaining -= edge;
            }
            else
            {
                l.C = c;
                l.HasEdge = true;
                remaining = 0;
            }
        }
    }

    /// <summary>Reports the left-maximal matches of this query start: every suffix whose common prefix with <c>query[start..]</c> reaches <paramref name="minLength"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private void Emit(in Locus m, int start, int[] qsym, int minLength, List<(int Position, int Length)> hits, List<(int PositionInText, int PositionInQuery, int Length)> results)
    {
        int lo = m.HasEdge ? m.C.Lb : m.V.Lb;
        int hi = m.HasEdge ? m.C.Rb : m.V.Rb;
        int len = m.Len;
        hits.Clear();
        for (int k = lo; k <= hi; k++) Consider(k, len, start, qsym, hits);

        int run = len;
        for (int k = lo; ; )
        {
            int v = _lcp[k];
            if (v < minLength) break;
            if (v < run) run = v;
            k--;
            Consider(k, run, start, qsym, hits);
        }

        run = len;
        for (int k = hi + 1; k <= _n; k++)
        {
            int v = _lcp[k];
            if (v < minLength) break;
            if (v < run) run = v;
            Consider(k, run, start, qsym, hits);
        }

        if (hits.Count > 1) hits.Sort(static (x, y) => x.Position.CompareTo(y.Position));
        foreach (var (position, length) in hits) results.Add((position, start, length));
        _emitted += hits.Count;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void Consider(int rank, int length, int start, int[] qsym, List<(int Position, int Length)> hits)
    {
        int r = _sa[rank];
        if (r == 0 || start == 0 || _text[r - 1] != qsym[start - 1]) hits.Add((r, length));
    }
}
