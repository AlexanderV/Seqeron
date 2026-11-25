using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>One alignment of a delta file (<c>DeltaAlignment_t</c> of <c>delta.hh</c>): reference / query coordinates (start &gt; end for the reverse strand), the error counts of the file, the percentages derived from them and the delta vector without its terminating zero... with it (the zero is stored as the last element, as the program's reader does).</summary>
public sealed class MummerDeltaAlignment
{
    /// <summary>Start in the reference.</summary>
    public long StartR { get; internal set; }

    /// <summary>End in the reference.</summary>
    public long EndR { get; internal set; }

    /// <summary>Start in the query.</summary>
    public long StartQ { get; internal set; }

    /// <summary>End in the query.</summary>
    public long EndQ { get; internal set; }

    /// <summary>Number of mismatches.</summary>
    public long IdentityErrors { get; internal set; }

    /// <summary>Number of similarity scores below 1.</summary>
    public long SimilarityErrors { get; internal set; }

    /// <summary>Number of stop codons / non-alphabetic characters.</summary>
    public long StopCount { get; internal set; }

    /// <summary>Percent identity (float, as the program computes it).</summary>
    public float Identity { get; internal set; }

    /// <summary>Percent similarity.</summary>
    public float Similarity { get; internal set; }

    /// <summary>Percent stop codons.</summary>
    public float Stop { get; internal set; }

    /// <summary>The delta vector including the terminating 0.</summary>
    public List<long> Deltas { get; } = [];
}

/// <summary>The alignments between one reference and one query record (<c>DeltaRecord_t</c>).</summary>
public sealed class MummerDeltaRecord
{
    /// <summary>Reference record id.</summary>
    public string ReferenceId { get; internal set; } = string.Empty;

    /// <summary>Query record id.</summary>
    public string QueryId { get; internal set; } = string.Empty;

    /// <summary>Reference length.</summary>
    public long ReferenceLength { get; internal set; }

    /// <summary>Query length.</summary>
    public long QueryLength { get; internal set; }

    /// <summary>The alignments.</summary>
    public List<MummerDeltaAlignment> Alignments { get; } = [];
}

/// <summary>A delta file (the output of <c>nucmer</c> / <c>promer</c>) read as MUMmer's <c>DeltaReader_t</c> reads it.</summary>
public sealed class MummerDeltaFile
{
    /// <summary>The reference path of the first line.</summary>
    public string ReferencePath { get; private set; } = string.Empty;

    /// <summary>The query path of the first line.</summary>
    public string QueryPath { get; private set; } = string.Empty;

    /// <summary><c>NUCMER</c> or <c>PROMER</c>.</summary>
    public string DataType { get; private set; } = string.Empty;

    /// <summary>The records in file order.</summary>
    public List<MummerDeltaRecord> Records { get; } = [];

    private sealed class Scanner(string text)
    {
        private int _pos;

        public int Peek() => _pos < text.Length ? text[_pos] : -1;

        public int Get() => _pos < text.Length ? text[_pos++] : -1;

        public string Token()
        {
            while (_pos < text.Length && IsSpace(text[_pos]))
                _pos++;
            int start = _pos;
            while (_pos < text.Length && !IsSpace(text[_pos]))
                _pos++;
            if (start == _pos)
                throw Fail();
            return text.Substring(start, _pos - start);
        }

        public long Long()
        {
            while (_pos < text.Length && IsSpace(text[_pos]))
                _pos++;
            int start = _pos;
            if (_pos < text.Length && (text[_pos] == '+' || text[_pos] == '-'))
                _pos++;
            int digits = _pos;
            while (_pos < text.Length && text[_pos] >= '0' && text[_pos] <= '9')
                _pos++;
            if (_pos == digits || !long.TryParse(text.AsSpan(start, _pos - start), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out long v))
                throw Fail();
            return v;
        }

        public bool AtEnd => _pos >= text.Length;

        private static bool IsSpace(char c) => c == ' ' || (c >= '\t' && c <= '\r');

        public static FormatException Fail() => new("ERROR: Could not parse delta file");
    }

    /// <summary>Parses the text of a delta file.</summary>
    /// <param name="text">The delta file.</param>
    /// <param name="readDeltas">False = <c>readNextHeadersOnly</c> (the delta vectors are not kept; <c>show-diff</c>).</param>
    /// <exception cref="FormatException">The program would print "ERROR: Could not parse delta file" and exit.</exception>
    public static MummerDeltaFile Parse(string text, bool readDeltas = true)
    {
        ArgumentNullException.ThrowIfNull(text);
        var sc = new Scanner(text);
        var file = new MummerDeltaFile
        {
            ReferencePath = sc.Token(),
            QueryPath = sc.Token(),
            DataType = sc.Token(),
        };
        if (file.DataType != "NUCMER" && file.DataType != "PROMER")
            throw Scanner.Fail();
        while (sc.Peek() != '>')
            if (sc.Get() == -1)
                break;

        while (sc.Peek() == '>')
        {
            sc.Get();
            var rec = new MummerDeltaRecord { ReferenceId = sc.Token(), QueryId = sc.Token(), ReferenceLength = sc.Long(), QueryLength = sc.Long() };
            if (rec.ReferenceLength <= 0 || rec.QueryLength <= 0)
                throw Scanner.Fail();
            while (sc.Get() != '\n')
                if (sc.AtEnd)
                    throw Scanner.Fail();
            while (sc.Peek() != '>' && sc.Peek() != -1)
            {
                var al = new MummerDeltaAlignment
                {
                    StartR = sc.Long(),
                    EndR = sc.Long(),
                    StartQ = sc.Long(),
                    EndQ = sc.Long(),
                    IdentityErrors = sc.Long(),
                    SimilarityErrors = sc.Long(),
                    StopCount = sc.Long(),
                };
                if (al.StartR <= 0 || al.EndR <= 0 || al.StartQ <= 0 || al.EndQ <= 0 || al.IdentityErrors < 0 || al.SimilarityErrors < 0 || al.StopCount < 0)
                    throw Scanner.Fail();
                float total = (float)(Math.Abs(al.EndR - al.StartR) + 1.0);
                if (file.DataType == "PROMER")
                    total = (float)(total / 3.0);
                long delta;
                do
                {
                    delta = sc.Long();
                    if (delta < 0)
                        total += 1;
                    if (readDeltas)
                        al.Deltas.Add(delta);
                }
                while (delta != 0);

                while (sc.Get() != '\n')
                    if (sc.AtEnd)
                        throw Scanner.Fail();
                float idyNum = total - al.IdentityErrors;
                float q = idyNum / total;
                al.Identity = (float)((double)q * 100.0);
                float simNum = total - al.SimilarityErrors;
                float qs = simNum / total;
                al.Similarity = (float)((double)qs * 100.0);
                al.Stop = (float)((double)(float)al.StopCount / ((double)total * 2.0) * 100.0);
                rec.Alignments.Add(al);
            }

            file.Records.Add(rec);
        }

        return file;
    }
}

/// <summary>A node of the <see cref="MummerDeltaGraph"/>: one reference or query record.</summary>
public sealed class MummerDeltaNode
{
    /// <summary>Record id.</summary>
    public string Id { get; internal set; } = string.Empty;

    /// <summary>Record length.</summary>
    public long Length { get; internal set; }

    /// <summary>The edges (one per delta record involving this node), in file order.</summary>
    public List<MummerDeltaEdge> Edges { get; } = [];

    /// <summary>The sequence, 1-based (index 0 and Length + 1 are zero), once loaded.</summary>
    public byte[]? Sequence { get; internal set; }

    internal int Order;
}

/// <summary>An SNP / indel found by <c>show-snps</c>.</summary>
public sealed class MummerSnp
{
    /// <summary>Distance to the nearest mismatch in the same alignment (<c>[BUFF]</c>).</summary>
    public long Buff { get; internal set; }

    /// <summary>The query / reference characters.</summary>
    public char CharQ { get; internal set; }

    /// <summary>The reference character.</summary>
    public char CharR { get; internal set; }

    /// <summary>The position in the query.</summary>
    public long PosQ { get; internal set; }

    /// <summary>The position in the reference.</summary>
    public long PosR { get; internal set; }

    /// <summary>Number of other alignments over the query position.</summary>
    public int ConflictQ { get; internal set; }

    /// <summary>Number of other alignments over the reference position.</summary>
    public int ConflictR { get; internal set; }

    /// <summary>Query context.</summary>
    public string ContextQ { get; internal set; } = string.Empty;

    /// <summary>Reference context.</summary>
    public string ContextR { get; internal set; } = string.Empty;

    internal MummerDeltaEdgelet Edgelet = null!;
    internal MummerDeltaEdge Edge = null!;
}

/// <summary>An alignment of the graph (<c>DeltaEdgelet_t</c>).</summary>
public sealed class MummerDeltaEdgelet
{
    /// <summary>The alignment meets the requirements (not removed by a filter).</summary>
    public bool IsGood { get; internal set; } = true;

    internal bool IsQLis, IsRLis, IsGLis;

    /// <summary>The reference strand is reverse.</summary>
    public bool ReverseR { get; internal set; }

    /// <summary>The query strand is reverse.</summary>
    public bool ReverseQ { get; internal set; }

    /// <summary>Identity (0 - 1).</summary>
    public float Identity { get; internal set; }

    /// <summary>Similarity (0 - 1).</summary>
    public float Similarity { get; internal set; }

    /// <summary>Stop codon fraction (0 - 1).</summary>
    public float Stop { get; internal set; }

    /// <summary>Counts of the file.</summary>
    public long IdentityErrors { get; internal set; }

    /// <summary>Similarity error count.</summary>
    public long SimilarityErrors { get; internal set; }

    /// <summary>Stop count (show-diff reuses it as a running index).</summary>
    public long StopCount { get; internal set; }

    /// <summary>Lower query bound.</summary>
    public long LoQ { get; internal set; }

    /// <summary>Upper query bound.</summary>
    public long HiQ { get; internal set; }

    /// <summary>Lower reference bound.</summary>
    public long LoR { get; internal set; }

    /// <summary>Upper reference bound.</summary>
    public long HiR { get; internal set; }

    /// <summary>Reading frames (set by show-snps).</summary>
    public int FrameQ { get; internal set; } = 1;

    /// <summary>Reading frame of the reference.</summary>
    public int FrameR { get; internal set; } = 1;

    /// <summary>The delta text, one value per line.</summary>
    public string Delta { get; internal set; } = string.Empty;

    /// <summary>The edge it belongs to.</summary>
    public MummerDeltaEdge? Edge { get; internal set; }

    internal List<MummerSnp> Snps = [];

    internal bool IsNegative => ReverseR != ReverseQ;

    internal bool IsPositive => ReverseR == ReverseQ;

    internal int Slope => ReverseR == ReverseQ ? 1 : -1;
}

/// <summary>The alignments of one reference / query pair.</summary>
public sealed class MummerDeltaEdge
{
    /// <summary>Reference node.</summary>
    public MummerDeltaNode ReferenceNode { get; internal set; } = null!;

    /// <summary>Query node.</summary>
    public MummerDeltaNode QueryNode { get; internal set; } = null!;

    /// <summary>The alignments.</summary>
    public List<MummerDeltaEdgelet> Edgelets { get; } = [];
}

/// <summary>
/// The delta graph of MUMmer 3.23 (<c>DeltaGraph_t</c>, <c>src/tigr/delta.cc</c>): reference and query records as nodes (ordered by id), a delta record as an edge, an alignment as an edgelet, and the filters used by <c>delta-filter</c>,
/// <c>show-diff</c>, <c>show-snps</c>: <c>flagScore</c> (-i / -l), <c>flagUNIQ</c> (-u), <c>flagQLIS</c> / <c>flagRLIS</c> (-q / -r, the longest increasing subset of alignments, the choice among equal ones made with glibc <c>rand()</c>, <c>srand(1)</c> as in the programs),
/// <c>flagGLIS</c> (-g), <c>flagMtoM</c> (-m), <c>flag1to1</c> (-1) and <c>clean</c> (with the unstable <c>std::partition</c> of libstdc++), then <c>outputDelta</c>.
/// </summary>
public sealed class MummerDeltaGraph
{
    private sealed class Lis
    {
        public MummerDeltaEdgelet A = null!;
        public long Score;
        public long Diff;
        public long From;
        public bool Used;
    }

    /// <summary>Reference nodes by id.</summary>
    public SortedDictionary<string, MummerDeltaNode> ReferenceNodes { get; } = new(StringComparer.Ordinal);

    /// <summary>Query nodes by id.</summary>
    public SortedDictionary<string, MummerDeltaNode> QueryNodes { get; } = new(StringComparer.Ordinal);

    /// <summary>Reference path of the delta file.</summary>
    public string ReferencePath { get; private set; } = string.Empty;

    /// <summary>Query path of the delta file.</summary>
    public string QueryPath { get; private set; } = string.Empty;

    /// <summary><c>NUCMER</c> or <c>PROMER</c>.</summary>
    public string DataType { get; private set; } = string.Empty;

    private readonly GlibcRandom _rand = new(1);

    /// <summary><c>DeltaGraph_t::build</c>.</summary>
    public static MummerDeltaGraph Build(string deltaText, bool readDeltas = true)
    {
        var file = MummerDeltaFile.Parse(deltaText, readDeltas);
        var g = new MummerDeltaGraph { ReferencePath = file.ReferencePath, QueryPath = file.QueryPath, DataType = file.DataType };
        int order = 0;
        foreach (var rec in file.Records)
        {
            var edge = new MummerDeltaEdge();
            if (!g.ReferenceNodes.TryGetValue(rec.ReferenceId, out var rn))
            {
                rn = new MummerDeltaNode { Id = rec.ReferenceId, Length = rec.ReferenceLength, Order = order++ };
                g.ReferenceNodes[rec.ReferenceId] = rn;
            }

            if (!g.QueryNodes.TryGetValue(rec.QueryId, out var qn))
            {
                qn = new MummerDeltaNode { Id = rec.QueryId, Length = rec.QueryLength, Order = order++ };
                g.QueryNodes[rec.QueryId] = qn;
            }

            edge.ReferenceNode = rn;
            edge.QueryNode = qn;
            foreach (var a in rec.Alignments)
            {
                var p = new MummerDeltaEdgelet
                {
                    Edge = edge,
                    Identity = (float)(a.Identity / 100.0),
                    Similarity = (float)(a.Similarity / 100.0),
                    Stop = (float)(a.Stop / 100.0),
                    IdentityErrors = a.IdentityErrors,
                    SimilarityErrors = a.SimilarityErrors,
                    StopCount = a.StopCount,
                    LoR = a.StartR,
                    HiR = a.EndR,
                    LoQ = a.StartQ,
                    HiQ = a.EndQ,
                };
                p.ReverseR = p.HiR < p.LoR;
                p.ReverseQ = p.HiQ < p.LoQ;
                var sb = new StringBuilder();
                foreach (long d in a.Deltas)
                    sb.Append(d.ToString(CultureInfo.InvariantCulture)).Append('\n');
                p.Delta = sb.ToString();
                if (p.ReverseR)
                    (p.LoR, p.HiR) = (p.HiR, p.LoR);
                if (p.ReverseQ)
                    (p.LoQ, p.HiQ) = (p.HiQ, p.LoQ);
                edge.Edgelets.Add(p);
            }

            rn.Edges.Add(edge);
            qn.Edges.Add(edge);
        }

        return g;
    }

    // std::partition (bidirectional version of libstdc++) with the predicate "not null": moves the survivors to the front, the order of the survivors is changed by the swaps
    private static int Partition<T>(IList<T?> a, Func<T?, bool> pred)
    {
        int first = 0, last = a.Count;
        while (true)
        {
            while (true)
            {
                if (first == last)
                    return first;
                if (pred(a[first]))
                    first++;
                else
                    break;
            }

            last--;
            while (true)
            {
                if (first == last)
                    return first;
                if (!pred(a[last]))
                    last--;
                else
                    break;
            }

            (a[first], a[last]) = (a[last], a[first]);
            first++;
        }
    }

    /// <summary><c>DeltaGraph_t::clean</c>: removes the alignments that are not good, then the edges and nodes left empty.</summary>
    public void Clean()
    {
        foreach (var key in ReferenceNodes.Keys.ToList())
        {
            var node = ReferenceNodes[key];
            var edges = new List<MummerDeltaEdge?>(node.Edges);
            foreach (var edge in edges.Select(e => e!).ToList())
            {
                var lets = new List<MummerDeltaEdgelet?>(edge!.Edgelets);
                for (int k = 0; k < lets.Count; k++)
                    if (!lets[k]!.IsGood)
                        lets[k] = null;
                int keep = Partition(lets, x => x != null);
                lets.RemoveRange(keep, lets.Count - keep);
                edge.Edgelets.Clear();
                edge.Edgelets.AddRange(lets!);
            }

            for (int j = 0; j < edges.Count; j++)
                if (edges[j]!.Edgelets.Count == 0)
                    edges[j] = null;
            int keepEdges = Partition(edges, x => x != null);
            edges.RemoveRange(keepEdges, edges.Count - keepEdges);
            node.Edges.Clear();
            node.Edges.AddRange(edges!);
            if (node.Edges.Count == 0)
                ReferenceNodes.Remove(key);
        }

        foreach (var key in QueryNodes.Keys.ToList())
        {
            var node = QueryNodes[key];
            var edges = new List<MummerDeltaEdge?>(node.Edges);
            for (int j = 0; j < edges.Count; j++)
                if (edges[j]!.Edgelets.Count == 0)
                    edges[j] = null;
            int keepEdges = Partition(edges, x => x != null);
            edges.RemoveRange(keepEdges, edges.Count - keepEdges);
            node.Edges.Clear();
            node.Edges.AddRange(edges!);
            if (node.Edges.Count == 0)
                QueryNodes.Remove(key);
        }
    }

    private static long ScoreLocal(long scoreJ, long leni, long lenj, long olap, float idyi, float maxolap)
    {
        if (olap > 0 && ((double)((float)olap / (float)leni) * 100.0 > maxolap || (double)((float)olap / (float)lenj) * 100.0 > maxolap))
            return -1;
        return scoreJ + (long)((leni - olap) * Math.Pow(idyi, 2));
    }

    private static long ScoreGlobal(long scoreJ, long leni, long olap, float idyi)
        => scoreJ + (long)((leni - olap) * Math.Pow(idyi, 2));

    private static long DiffAligns(MummerDeltaEdgelet i, MummerDeltaEdgelet j)
        => (j.LoR < i.LoR ? Math.Abs(j.HiR - i.LoR) : Math.Abs(i.HiR - j.LoR)) + (j.LoQ < i.LoQ ? Math.Abs(j.HiQ - i.LoQ) : Math.Abs(i.HiQ - j.LoQ));

    private static bool QLess(MummerDeltaEdgelet i, MummerDeltaEdgelet j)
    {
        if (i.LoQ < j.LoQ)
            return true;
        if (i.LoQ > j.LoQ)
            return false;
        return ScoreLocal(0, i.HiQ - i.LoQ + 1, 0, 0, i.Identity, 0) > ScoreLocal(0, j.HiQ - j.LoQ + 1, 0, 0, j.Identity, 0);
    }

    private static bool RLess(MummerDeltaEdgelet i, MummerDeltaEdgelet j)
    {
        if (i.LoR < j.LoR)
            return true;
        if (i.LoR > j.LoR)
            return false;
        return ScoreLocal(0, i.HiR - i.LoR + 1, 0, 0, i.Identity, 0) > ScoreLocal(0, j.HiR - j.LoR + 1, 0, 0, j.Identity, 0);
    }

    private long PickBest(Lis[] lis, List<long> allbest, float epsilon)
    {
        long size = allbest.Count;
        if (epsilon < 0 && size != 0)
        {
            long eqc = 0;
            for (; eqc < size; ++eqc)
                if (lis[allbest[(int)eqc]].Diff != lis[allbest[0]].Diff)
                    break;
            return (int)((double)eqc * _rand.Next() / (2147483647 + 1.0));
        }

        return size;
    }

    private static bool UpdateBest(Lis[] lis, long size, List<long> allbest, float epsilon)
    {
        if (size == 0)
            return false;
        long best, i;
        for (best = 0; best < size; ++best)
            if (!lis[best].Used)
                break;
        for (i = best + 1; i < size; ++i)
            if (!lis[i].Used && (lis[i].Score > lis[best].Score || (lis[i].Score == lis[best].Score && lis[i].Diff < lis[best].Diff)))
                best = i;

        if (allbest.Count != 0
            && (best == size
                || (epsilon < 0 && lis[best].Score < lis[allbest[0]].Score)
                || (epsilon >= 0 && (double)((float)(lis[allbest[0]].Score - lis[best].Score) / (float)lis[allbest[0]].Score) * 100.0 > epsilon)))
            return false;

        allbest.Add(best);
        for (i = best; i >= 0 && i < size; i = lis[i].From)
            lis[i].Used = true;
        return true;
    }

    /// <summary><c>flagGOOD</c>.</summary>
    public void FlagGood()
    {
        foreach (var n in ReferenceNodes.Values)
            foreach (var e in n.Edges)
                foreach (var l in e.Edgelets)
                    l.IsGood = true;
    }

    /// <summary><c>flag1to1</c>.</summary>
    public void Flag1To1(float epsilon = -1, float maxOverlap = 100f)
    {
        FlagRLis(epsilon, maxOverlap, false);
        FlagQLis(epsilon, maxOverlap, false);
        foreach (var n in ReferenceNodes.Values)
            foreach (var e in n.Edges)
                foreach (var l in e.Edgelets)
                    if (!l.IsRLis || !l.IsQLis)
                        l.IsGood = false;
    }

    /// <summary><c>flagMtoM</c>.</summary>
    public void FlagMtoM(float epsilon = -1, float maxOverlap = 100f)
    {
        FlagRLis(epsilon, maxOverlap, false);
        FlagQLis(epsilon, maxOverlap, false);
        foreach (var n in ReferenceNodes.Values)
            foreach (var e in n.Edges)
                foreach (var l in e.Edgelets)
                    if (!l.IsRLis && !l.IsQLis)
                        l.IsGood = false;
    }

    /// <summary><c>flagGLIS</c> (<c>delta-filter -g</c>): the longest global chain of mutually consistent alignments per reference / query pair.</summary>
    public void FlagGlis(float epsilon = -1)
    {
        foreach (var mn in ReferenceNodes.Values)
        {
            foreach (var ei in mn.Edges)
            {
                var edgelets = new List<MummerDeltaEdgelet>();
                foreach (var el in ei.Edgelets)
                {
                    if (!el.IsGood)
                        continue;
                    edgelets.Add(el);
                    if (el.ReverseR == el.ReverseQ)
                    {
                        el.ReverseQ = false;
                    }
                    else
                    {
                        if (el.ReverseQ)
                            (el.LoQ, el.HiQ) = (el.HiQ, el.LoQ);
                        el.LoQ = ei.QueryNode.Length - el.LoQ + 1;
                        el.HiQ = ei.QueryNode.Length - el.HiQ + 1;
                        el.ReverseQ = true;
                    }
                }

                int n = edgelets.Count;
                var lis = new Lis[n];
                for (int k = 0; k < n; k++)
                    lis[k] = new Lis();
                LibstdcxxSort.Sort(edgelets, QLess);
                var allbest = new List<long>();
                do
                {
                    for (int i = 0; i < n; ++i)
                    {
                        if (lis[i].Used)
                            continue;
                        lis[i].A = edgelets[i];
                        long lenR = lis[i].A.HiR - lis[i].A.LoR + 1;
                        long lenQ = lis[i].A.HiQ - lis[i].A.LoQ + 1;
                        long len = lenR > lenQ ? lenQ : lenR;
                        lis[i].Score = ScoreGlobal(0, len, 0, lis[i].A.Identity);
                        lis[i].From = -1;
                        lis[i].Diff = 0;
                        for (int j = 0; j < i; ++j)
                        {
                            if (lis[j].Used)
                                continue;
                            if (lis[i].A.ReverseQ != lis[j].A.ReverseQ)
                                continue;
                            lenR = lis[i].A.HiR - lis[i].A.LoR + 1;
                            lenQ = lis[i].A.HiQ - lis[i].A.LoQ + 1;
                            len = lenR > lenQ ? lenQ : lenR;
                            long olapR = lis[j].A.HiR - lis[i].A.LoR + 1;
                            long olapQ = lis[j].A.HiQ - lis[i].A.LoQ + 1;
                            long olap = olapR > olapQ ? olapR : olapQ;
                            if (olap < 0)
                                olap = 0;
                            long diff = lis[j].Diff + DiffAligns(lis[i].A, lis[j].A);
                            long score = ScoreGlobal(lis[j].Score, len, olap, lis[i].A.Identity);
                            if (score > lis[i].Score || (score == lis[i].Score && diff < lis[i].Diff))
                            {
                                lis[i].From = j;
                                lis[i].Score = score;
                                lis[i].Diff = diff;
                            }
                        }
                    }
                }
                while (UpdateBest(lis, n, allbest, epsilon));

                MarkChains(lis, n, allbest, epsilon, (a) => a.IsGLis = true);

                foreach (var el in edgelets)
                {
                    if (!el.IsGLis)
                        el.IsGood = false;
                    if (!el.ReverseQ)
                    {
                        el.ReverseQ = el.ReverseR;
                    }
                    else
                    {
                        if (!el.ReverseR)
                            (el.LoQ, el.HiQ) = (el.HiQ, el.LoQ);
                        el.LoQ = ei.QueryNode.Length - el.LoQ + 1;
                        el.HiQ = ei.QueryNode.Length - el.HiQ + 1;
                        el.ReverseQ = !el.ReverseR;
                    }
                }
            }
        }
    }

    private void MarkChains(Lis[] lis, int n, List<long> allbest, float epsilon, Action<MummerDeltaEdgelet> mark)
    {
        long beg = PickBest(lis, allbest, epsilon);
        long end = allbest.Count;
        if (beg == end)
            beg = 0;
        else
            end = beg + 1;
        for (; beg < end; ++beg)
            for (long i = allbest[(int)beg]; i >= 0 && i < n; i = lis[i].From)
                mark(lis[i].A);
    }

    /// <summary><c>flagQLIS</c> (<c>delta-filter -q</c>).</summary>
    public void FlagQLis(float epsilon, float maxOverlap, bool flagBad = true)
        => FlagLis(QueryNodes, epsilon, maxOverlap, flagBad, true);

    /// <summary><c>flagRLIS</c> (<c>delta-filter -r</c>).</summary>
    public void FlagRLis(float epsilon, float maxOverlap, bool flagBad = true)
        => FlagLis(ReferenceNodes, epsilon, maxOverlap, flagBad, false);

    private void FlagLis(SortedDictionary<string, MummerDeltaNode> nodes, float epsilon, float maxOverlap, bool flagBad, bool query)
    {
        foreach (var mn in nodes.Values)
        {
            var edgelets = new List<MummerDeltaEdgelet>();
            foreach (var ei in mn.Edges)
                foreach (var el in ei.Edgelets)
                    if (el.IsGood)
                        edgelets.Add(el);

            int n = edgelets.Count;
            var lis = new Lis[n];
            for (int k = 0; k < n; k++)
                lis[k] = new Lis();
            LibstdcxxSort.Sort(edgelets, query ? QLess : RLess);

            long Lo(MummerDeltaEdgelet a) => query ? a.LoQ : a.LoR;
            long Hi(MummerDeltaEdgelet a) => query ? a.HiQ : a.HiR;

            var allbest = new List<long>();
            do
            {
                for (int i = 0; i < n; ++i)
                {
                    if (lis[i].Used)
                        continue;
                    lis[i].A = edgelets[i];
                    long leni = Hi(lis[i].A) - Lo(lis[i].A) + 1;
                    lis[i].Score = ScoreLocal(0, leni, 0, 0, lis[i].A.Identity, 0);
                    lis[i].From = -1;
                    lis[i].Diff = 0;
                    for (int j = 0; j < i; ++j)
                    {
                        if (lis[j].Used)
                            continue;
                        if (lis[j].From >= 0 && Hi(lis[lis[j].From].A) >= Lo(lis[i].A))
                            continue;
                        leni = Hi(lis[i].A) - Lo(lis[i].A) + 1;
                        long lenj = Hi(lis[j].A) - Lo(lis[j].A) + 1;
                        long olap = Hi(lis[j].A) - Lo(lis[i].A) + 1;
                        if (olap < 0)
                            olap = 0;
                        long diff = lis[j].Diff + DiffAligns(lis[i].A, lis[j].A);
                        long score = ScoreLocal(lis[j].Score, leni, lenj, olap, lis[i].A.Identity, maxOverlap);
                        if (score > lis[i].Score || (score == lis[i].Score && diff < lis[i].Diff))
                        {
                            lis[i].From = j;
                            lis[i].Score = score;
                            lis[i].Diff = diff;
                        }
                    }
                }
            }
            while (UpdateBest(lis, n, allbest, epsilon));

            if (query)
                MarkChains(lis, n, allbest, epsilon, a => a.IsQLis = true);
            else
                MarkChains(lis, n, allbest, epsilon, a => a.IsRLis = true);

            if (flagBad)
                foreach (var el in edgelets)
                    if (query ? !el.IsQLis : !el.IsRLis)
                        el.IsGood = false;
        }
    }

    /// <summary><c>flagScore</c> (<c>delta-filter -i -l</c>).</summary>
    public void FlagScore(long minLength, float minIdentity)
    {
        foreach (var n in ReferenceNodes.Values)
            foreach (var e in n.Edges)
                foreach (var l in e.Edgelets)
                {
                    if (!l.IsGood)
                        continue;
                    if ((double)l.Identity * 100.0 < minIdentity)
                        l.IsGood = false;
                    if (l.HiR - l.LoR + 1 < minLength || l.HiQ - l.LoQ + 1 < minLength)
                        l.IsGood = false;
                }
    }

    /// <summary><c>flagUNIQ</c> (<c>delta-filter -u</c>).</summary>
    public void FlagUniq(float minUnique)
    {
        foreach (bool query in new[] { false, true })
        {
            foreach (var mn in (query ? QueryNodes : ReferenceNodes).Values)
            {
                long len0 = mn.Length;
                var cov = new byte[len0 + 2];
                var edgelets = new List<MummerDeltaEdgelet>();
                foreach (var ei in mn.Edges)
                    foreach (var el in ei.Edgelets)
                        if (el.IsGood)
                        {
                            edgelets.Add(el);
                            long lo = query ? el.LoQ : el.LoR, hi = query ? el.HiQ : el.HiR;
                            for (long i = lo; i <= hi; i++)
                                if (cov[i] < byte.MaxValue)
                                    cov[i]++;
                        }

                foreach (var el in edgelets)
                {
                    long uniq = 0;
                    long lo = query ? el.LoQ : el.LoR, hi = query ? el.HiQ : el.HiR;
                    long len = hi - lo + 1;
                    for (long i = lo; i <= hi; i++)
                        if (cov[i] == 1)
                            uniq++;
                    if ((double)((float)uniq / (float)len) * 100.0 < minUnique)
                        el.IsGood = false;
                }
            }
        }
    }

    /// <summary><c>loadSequences</c>: attaches the sequences of the FASTA files named in the delta header to the nodes.</summary>
    /// <exception cref="FormatException">A record of the delta file is missing from a FASTA file or has another length (the program's error exits).</exception>
    public void LoadSequences(string referenceFasta, string queryFasta)
    {
        Load(referenceFasta, ReferenceNodes, "Reference");
        Load(queryFasta, QueryNodes, "Query");
        foreach (var n in ReferenceNodes)
            if (n.Value.Sequence == null)
                throw new FormatException("ERROR: '" + n.Key + "' not found in reference file");
        foreach (var n in QueryNodes)
            if (n.Value.Sequence == null)
                throw new FormatException("ERROR: '" + n.Key + "' not found in query file");

        static void Load(string fasta, SortedDictionary<string, MummerDeltaNode> nodes, string what)
        {
            foreach (var (id, seq) in NucmerPostnuc.ReadFasta(fasta, out _))
            {
                if (!nodes.TryGetValue(id, out var node))
                    continue;
                var bytes = new byte[seq.Length + 2];
                for (int i = 0; i < seq.Length; i++)
                    bytes[i + 1] = (byte)seq[i];
                node.Sequence = bytes;
                if (seq.Length != node.Length)
                    throw new FormatException("ERROR: " + what + " input does not match delta file");
            }
        }
    }

    /// <summary><c>outputDelta</c>: the delta text of the good alignments, grouped by query node.</summary>
    public string OutputDelta()
    {
        var sb = new StringBuilder();
        sb.Append(ReferencePath).Append(' ').Append(QueryPath).Append('\n').Append(DataType == "PROMER" ? "PROMER" : "NUCMER").Append('\n');
        foreach (var mn in QueryNodes.Values)
        {
            foreach (var ei in mn.Edges)
            {
                bool header = false;
                foreach (var el in ei.Edgelets)
                {
                    if (!el.IsGood)
                        continue;
                    if (!header)
                    {
                        sb.Append('>').Append(ei.ReferenceNode.Id).Append(' ').Append(ei.QueryNode.Id).Append(' ')
                            .Append(ei.ReferenceNode.Length.ToString(CultureInfo.InvariantCulture)).Append(' ').Append(ei.QueryNode.Length.ToString(CultureInfo.InvariantCulture)).Append('\n');
                        header = true;
                    }

                    long s1 = el.LoR, e1 = el.HiR, s2 = el.LoQ, e2 = el.HiQ;
                    if (el.ReverseR)
                        (s1, e1) = (e1, s1);
                    if (el.ReverseQ)
                        (s2, e2) = (e2, s2);
                    var ci = CultureInfo.InvariantCulture;
                    sb.Append(s1.ToString(ci)).Append(' ').Append(e1.ToString(ci)).Append(' ').Append(s2.ToString(ci)).Append(' ').Append(e2.ToString(ci)).Append(' ')
                        .Append(el.IdentityErrors.ToString(ci)).Append(' ').Append(el.SimilarityErrors.ToString(ci)).Append(' ').Append(el.StopCount.ToString(ci)).Append('\n').Append(el.Delta);
                }
            }
        }

        return sb.ToString();
    }
}
