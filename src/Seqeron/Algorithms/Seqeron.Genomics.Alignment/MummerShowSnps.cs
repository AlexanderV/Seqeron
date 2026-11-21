using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>Options of <see cref="MummerShowSnps.Run"/> (<c>show-snps</c> 3.23).</summary>
public sealed record MummerShowSnpsOptions
{
    /// <summary><c>-C</c> inverted: report SNPs of ambiguous alignments and show the [R] [Q] columns (default true).</summary>
    public bool ShowConflict { get; init; } = true;

    /// <summary><c>-I</c> inverted: report indels (default true).</summary>
    public bool ShowIndels { get; init; } = true;

    /// <summary><c>-l</c>: sequence length columns.</summary>
    public bool ShowLength { get; init; }

    /// <summary><c>-q</c>: sort by query (ignored with -r).</summary>
    public bool SortQuery { get; init; }

    /// <summary><c>-r</c>: sort by reference (the default when neither is given).</summary>
    public bool SortReference { get; init; }

    /// <summary><c>-T</c>: tab-delimited.</summary>
    public bool Tabular { get; init; }

    /// <summary><c>-H</c> inverted: print the header (default true).</summary>
    public bool PrintHeader { get; init; } = true;

    /// <summary><c>-x</c>: characters of SNP context on each side.</summary>
    public int Context { get; init; }

    /// <summary><c>-S</c>: only the alignments of these <c>show-coords</c> lines (the program reads them from standard input); null = all.</summary>
    public IReadOnlyList<string>? SelectedAlignments { get; init; }
}

/// <summary>
/// <c>show-snps</c> of MUMmer 3.23 (<c>src/tigr/show-snps.cc</c>): the single-base substitutions and indels of every alignment of a nucmer delta file, with their positions in both sequences, the distance to the
/// nearest other difference of the alignment ([BUFF]), to the nearest sequence end ([DIST]), the number of other alignments over the position ([R] [Q]) and optional sequence context.
/// </summary>
public static class MummerShowSnps
{
    private const char IndelChar = '.';
    private const char SeqEndChar = '-';

    private static bool Iupac(char xc, char yc)
    {
        char x = char.ToLowerInvariant(xc), y = char.ToLowerInvariant(yc);
        if (x == 'n' || x == 'x' || y == 'n' || y == 'x')
            return true;
        string set = x switch
        {
            'a' => "arwmdhv",
            'c' => "cysmbhv",
            'g' => "grskbdv",
            't' => "tywkbdh",
            'r' => "ragdv",
            'y' => "yctbh",
            's' => "scgbv",
            'w' => "watdh",
            'm' => "machv",
            'k' => "kgtbd",
            'b' => "bcgt",
            'd' => "dagt",
            'h' => "hact",
            'v' => "vacg",
            _ => string.Empty,
        };
        return set.Contains(y, StringComparison.Ordinal);
    }

    private static long RevC(long coord, long len) => len - coord + 1;

    private static long Norm(long c, long l, int f) => f < 0 ? RevC(c, l) : c;

    private static bool SnpRLess(MummerSnp a, MummerSnp b)
    {
        int i = string.CompareOrdinal(a.Edge.ReferenceNode.Id, b.Edge.ReferenceNode.Id);
        if (i < 0)
            return true;
        if (i > 0)
            return false;
        if (a.PosR < b.PosR)
            return true;
        if (a.PosR > b.PosR)
            return false;
        int j = string.CompareOrdinal(a.Edge.QueryNode.Id, b.Edge.QueryNode.Id);
        if (j < 0)
            return true;
        if (j > 0)
            return false;
        return a.PosQ < b.PosQ;
    }

    private static bool SnpQLess(MummerSnp a, MummerSnp b)
    {
        int i = string.CompareOrdinal(a.Edge.QueryNode.Id, b.Edge.QueryNode.Id);
        if (i < 0)
            return true;
        if (i > 0)
            return false;
        if (a.PosQ < b.PosQ)
            return true;
        if (a.PosQ > b.PosQ)
            return false;
        int j = string.CompareOrdinal(a.Edge.ReferenceNode.Id, b.Edge.ReferenceNode.Id);
        if (j < 0)
            return true;
        if (j > 0)
            return false;
        return a.PosR < b.PosR;
    }

    /// <summary>Runs <c>show-snps</c>.</summary>
    /// <param name="deltaText">The nucmer delta file.</param>
    /// <param name="referenceFasta">The reference FASTA named in the delta header.</param>
    /// <param name="queryFasta">The query FASTA named in the delta header.</param>
    /// <param name="options">Switches (null = defaults).</param>
    /// <exception cref="FormatException">The delta file cannot be parsed, a record is missing from a FASTA file, or a -S alignment is not found.</exception>
    /// <exception cref="NotSupportedException">The delta file is a promer file (a registered sub-item).</exception>
    public static MummerTextOutput Run(string deltaText, string referenceFasta, string queryFasta, MummerShowSnpsOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(deltaText);
        ArgumentNullException.ThrowIfNull(referenceFasta);
        ArgumentNullException.ThrowIfNull(queryFasta);
        var o = options ?? new MummerShowSnpsOptions();
        var err = new StringBuilder();
        if (o.Context < 0)
            throw new ArgumentOutOfRangeException(nameof(options), "ERROR: SNP context must be a positive int");
        bool sortReference = o.SortReference;
        if (o.SortReference && o.SortQuery)
            err.Append("WARNING: both -r and -q were passed, -q ignored\n");
        if (!o.SortReference && !o.SortQuery)
            sortReference = true;

        HashSet<string>? select = null;
        if (o.SelectedAlignments != null)
            select = SelectAligns(o.SelectedAlignments);

        var graph = MummerDeltaGraph.Build(deltaText, true);
        if (graph.DataType != "NUCMER")
            throw new NotSupportedException("show-snps of a promer delta file (amino-acid substitutions) is not ported.");
        graph.LoadSequences(referenceFasta, queryFasta);
        FindSnps(graph, o, sortReference, select);
        if (select != null && select.Count != 0)
            throw new FormatException("ERROR: One or more alignments from stdin could not be found");
        CheckSnps(graph);

        var snps = new List<MummerSnp>();
        foreach (var mn in graph.ReferenceNodes.Values)
            foreach (var ei in mn.Edges)
                foreach (var li in ei.Edgelets)
                    foreach (var si in li.Snps)
                        if ((o.ShowConflict || (si.ConflictR == 0 && si.ConflictQ == 0)) && (o.ShowIndels || (si.CharR != IndelChar && si.CharQ != IndelChar)))
                            snps.Add(si);

        if (sortReference)
            LibstdcxxSort.Sort(snps, SnpRLess);
        else
            LibstdcxxSort.Sort(snps, SnpQLess);

        var sb = new StringBuilder();
        if (o.Tabular)
            PrintTabular(sb, snps, graph, o);
        else
            PrintHuman(sb, snps, graph, o);
        return new MummerTextOutput(sb.ToString(), err.ToString());
    }

    private static HashSet<string> SelectAligns(IReadOnlyList<string> lines)
    {
        var set = new HashSet<string>(StringComparer.Ordinal);
        foreach (string line in lines)
        {
            if (line.Length == 0)
                continue;
            var tok = line.Split([' ', '\t', '\n', '\r', '\v', '\f'], StringSplitOptions.RemoveEmptyEntries).ToList();
            int p = 0;
            string s1, e1, s2, e2, tR, tQ;
            if (tok.Count < 3)
                throw new FormatException("ERROR: Could not parse input from stdin");
            s1 = tok[p++];
            e1 = tok[p++];
            s2 = tok[p++];
            if (s2 == "|")
            {
                if (p >= tok.Count)
                    throw new FormatException("ERROR: Could not parse input from stdin");
                s2 = tok[p++];
            }

            if (p + 3 > tok.Count)
                throw new FormatException("ERROR: Could not parse input from stdin");
            e2 = tok[p++];
            tR = tok[p++];
            tQ = tok[p++];
            while (p < tok.Count)
            {
                tR = tQ;
                tQ = tok[p++];
            }

            set.Add(string.Join(' ', s1, e1, s2, e2, tR, tQ));
        }

        return set;
    }

    private static void FindSnps(MummerDeltaGraph graph, MummerShowSnpsOptions o, bool sortReference, HashSet<string>? select)
    {
        var ci = CultureInfo.InvariantCulture;
        foreach (var mn in graph.ReferenceNodes.Values)
        {
            foreach (var ei in mn.Edges)
            {
                var R = new byte[7][];
                var Q = new byte[7][];
                R[0] = ei.ReferenceNode.Sequence!;
                Q[0] = ei.QueryNode.Sequence!;
                long lenR = ei.ReferenceNode.Length, lenQ = ei.QueryNode.Length;

                foreach (var li in ei.Edgelets)
                {
                    long alenR = lenR, alenQ = lenQ;
                    if (select != null)
                    {
                        string key = string.Concat(
                            (li.ReverseR ? li.HiR : li.LoR).ToString(ci), " ", (li.ReverseR ? li.LoR : li.HiR).ToString(ci), " ",
                            (li.ReverseQ ? li.HiQ : li.LoQ).ToString(ci), " ", (li.ReverseQ ? li.LoQ : li.HiQ).ToString(ci), " ",
                            ei.ReferenceNode.Id, " ", ei.QueryNode.Id);
                        if (!select.Remove(key))
                            continue;
                    }

                    long sR, eR, sQ;
                    int frameR = 1, frameQ = 1;
                    if (li.ReverseR)
                    {
                        sR = RevC(li.HiR, lenR);
                        eR = RevC(li.LoR, lenR);
                        frameR += 3;
                    }
                    else
                    {
                        sR = li.LoR;
                        eR = li.HiR;
                    }

                    if (li.ReverseQ)
                    {
                        sQ = RevC(li.HiQ, lenQ);
                        frameQ += 3;
                    }
                    else
                    {
                        sQ = li.LoQ;
                    }

                    int ri = frameR, qi = frameQ;
                    if (frameR > 3)
                        frameR = -(frameR - 3);
                    if (frameQ > 3)
                        frameQ = -(frameQ - 3);

                    if (R[ri] == null)
                        R[ri] = li.ReverseR ? MummerSequence.ReverseComplement(R[0]) : (byte[])R[0].Clone();
                    if (Q[qi] == null)
                        Q[qi] = li.ReverseQ ? MummerSequence.ReverseComplement(Q[0]) : (byte[])Q[0].Clone();

                    long rpos = sR, qpos = sQ, remain = eR - sR + 1;
                    li.FrameR = frameR;
                    li.FrameQ = frameQ;
                    var rs = R[ri];
                    var qs = Q[qi];

                    MummerSnp Substitution()
                    {
                        var snp = new MummerSnp { Edge = ei, Edgelet = li };
                        snp.PosR = Norm(rpos - 1, lenR, frameR);
                        snp.PosQ = Norm(qpos - 1, lenQ, frameQ);
                        snp.CharR = char.ToUpperInvariant((char)rs[rpos - 1]);
                        snp.CharQ = char.ToUpperInvariant((char)qs[qpos - 1]);
                        var cr = new StringBuilder();
                        for (long rctx = rpos - o.Context - 1; rctx < rpos + o.Context; rctx++)
                        {
                            if (rctx < 1 || rctx > alenR)
                                cr.Append(SeqEndChar);
                            else if (rctx == rpos - 1)
                                cr.Append(snp.CharR);
                            else
                                cr.Append(char.ToUpperInvariant((char)rs[rctx]));
                        }

                        var cq = new StringBuilder();
                        for (long qctx = qpos - o.Context - 1; qctx < qpos + o.Context; qctx++)
                        {
                            if (qctx < 1 || qctx > alenQ)
                                cq.Append(SeqEndChar);
                            else if (qctx == qpos - 1)
                                cq.Append(snp.CharQ);
                            else
                                cq.Append(char.ToUpperInvariant((char)qs[qctx]));
                        }

                        snp.ContextR = cr.ToString();
                        snp.ContextQ = cq.ToString();
                        return snp;
                    }

                    using var reader = new StringReader(li.Delta);
                    string? dl;
                    while ((dl = reader.ReadLine()) != null)
                    {
                        long delta = long.Parse(dl, ci);
                        if (delta == 0)
                            break;
                        int sign = delta > 0 ? 1 : -1;
                        delta = Math.Abs(delta);
                        long i;
                        for (i = 1; i < delta; i++)
                        {
                            if (rs[rpos++] != qs[qpos++])
                            {
                                if (Iupac((char)rs[rpos - 1], (char)qs[qpos - 1]))
                                    continue;
                                li.Snps.Add(Substitution());
                            }
                        }

                        var snp = new MummerSnp { Edge = ei, Edgelet = li };
                        var cr = new StringBuilder();
                        long rctx2, qctx2;
                        for (rctx2 = rpos - o.Context; rctx2 < rpos; rctx2++)
                            cr.Append(rctx2 < 1 ? SeqEndChar : char.ToUpperInvariant((char)rs[rctx2]));
                        var cq = new StringBuilder();
                        for (qctx2 = qpos - o.Context; qctx2 < qpos; qctx2++)
                            cq.Append(qctx2 < 1 ? SeqEndChar : char.ToUpperInvariant((char)qs[qctx2]));

                        if (sign > 0)
                        {
                            snp.PosR = Norm(rpos, lenR, frameR);
                            snp.PosQ = frameQ > 0 ? Norm(qpos - 1, lenQ, frameQ) : Norm(qpos, lenQ, frameQ);
                            snp.CharR = char.ToUpperInvariant((char)rs[rpos++]);
                            snp.CharQ = IndelChar;
                            remain -= i;
                            rctx2++;
                        }
                        else
                        {
                            snp.PosQ = Norm(qpos, lenQ, frameQ);
                            snp.PosR = frameR > 0 ? Norm(rpos - 1, lenR, frameR) : Norm(rpos, lenR, frameR);
                            snp.CharR = IndelChar;
                            snp.CharQ = char.ToUpperInvariant((char)qs[qpos++]);
                            remain -= i - 1;
                            qctx2++;
                        }

                        cr.Append(snp.CharR);
                        for (; rctx2 < rpos + o.Context; rctx2++)
                            cr.Append(rctx2 > alenR ? SeqEndChar : char.ToUpperInvariant((char)rs[rctx2]));
                        cq.Append(snp.CharQ);
                        for (; qctx2 < qpos + o.Context; qctx2++)
                            cq.Append(qctx2 > alenQ ? SeqEndChar : char.ToUpperInvariant((char)qs[qctx2]));
                        snp.ContextR = cr.ToString();
                        snp.ContextQ = cq.ToString();
                        li.Snps.Add(snp);
                    }

                    for (long i = 0; i < remain; i++)
                    {
                        if (rs[rpos++] != qs[qpos++])
                        {
                            if (Iupac((char)rs[rpos - 1], (char)qs[qpos - 1]))
                                continue;
                            li.Snps.Add(Substitution());
                        }
                    }

                    if (sortReference)
                    {
                        LibstdcxxSort.Sort(li.Snps, SnpRLess);
                        for (int si = 0; si < li.Snps.Count; si++)
                        {
                            var s = li.Snps[si];
                            s.Buff = 1 + (s.PosR - li.LoR < li.HiR - s.PosR ? s.PosR - li.LoR : li.HiR - s.PosR);
                            if (si - 1 >= 0 && s.PosR - li.Snps[si - 1].PosR < s.Buff)
                                s.Buff = s.PosR - li.Snps[si - 1].PosR;
                            if (si + 1 < li.Snps.Count && li.Snps[si + 1].PosR - s.PosR < s.Buff)
                                s.Buff = li.Snps[si + 1].PosR - s.PosR;
                        }
                    }
                    else
                    {
                        LibstdcxxSort.Sort(li.Snps, SnpQLess);
                        for (int si = 0; si < li.Snps.Count; si++)
                        {
                            var s = li.Snps[si];
                            s.Buff = 1 + (s.PosQ - li.LoQ < li.HiQ - s.PosQ ? s.PosQ - li.LoQ : li.HiQ - s.PosQ);
                            if (si - 1 >= 0 && s.PosQ - li.Snps[si - 1].PosQ < s.Buff)
                                s.Buff = s.PosQ - li.Snps[si - 1].PosQ;
                            if (si + 1 < li.Snps.Count && li.Snps[si + 1].PosQ - s.PosQ < s.Buff)
                                s.Buff = li.Snps[si + 1].PosQ - s.PosQ;
                        }
                    }
                }
            }
        }
    }

    private static void CheckSnps(MummerDeltaGraph graph)
    {
        foreach (bool query in new[] { false, true })
        {
            foreach (var mn in (query ? graph.QueryNodes : graph.ReferenceNodes).Values)
            {
                var cov = new byte[mn.Length + 2];
                foreach (var ei in mn.Edges)
                    foreach (var el in ei.Edgelets)
                        for (long i = query ? el.LoQ : el.LoR; i <= (query ? el.HiQ : el.HiR); i++)
                            if (cov[i] < byte.MaxValue)
                                cov[i]++;
                foreach (var ei in mn.Edges)
                    foreach (var el in ei.Edgelets)
                        foreach (var si in el.Snps)
                        {
                            if (query)
                                si.ConflictQ = cov[si.PosQ] - 1;
                            else
                                si.ConflictR = cov[si.PosR] - 1;
                        }
            }
        }
    }

    private static string L(long v, int w = 0) => v.ToString(CultureInfo.InvariantCulture).PadLeft(w);

    private static void PrintHuman(StringBuilder sb, List<MummerSnp> snps, MummerDeltaGraph graph, MummerShowSnpsOptions o)
    {
        int ctxw = (2 * o.Context) + 1;
        int ctxc = ctxw < 7 ? 7 : ctxw;
        if (o.PrintHeader)
        {
            sb.Append(graph.ReferencePath).Append(' ').Append(graph.QueryPath).Append('\n').Append(graph.DataType).Append("\n\n");
            sb.Append("[P1]".PadLeft(8)).Append("  ").Append("[SUB]".PadLeft(5)).Append("  ").Append("[P2]".PadRight(8)).Append("  | ");
            sb.Append("[BUFF]".PadLeft(8)).Append(' ').Append("[DIST]".PadLeft(8)).Append("  | ");
            if (o.ShowConflict)
                sb.Append("[R]".PadLeft(4)).Append(' ').Append("[Q]".PadLeft(4)).Append("  | ");
            if (o.ShowLength)
                sb.Append("[LEN R]".PadLeft(8)).Append(' ').Append("[LEN Q]".PadLeft(8)).Append("  | ");
            if (o.Context != 0)
            {
                sb.Append(' ', ctxc - 7).Append(" [CTX R]  ").Append(' ', ctxc - 7).Append("[CTX Q]  | ");
            }

            sb.Append("[FRM]".PadLeft(5)).Append("  ").Append("[TAGS]").Append('\n');
            if (o.ShowConflict)
                sb.Append("=============");
            if (o.ShowLength)
                sb.Append("=====================");
            if (o.Context != 0)
                sb.Append('=', (2 * ctxc) + 7);
            sb.Append("=================================").Append("==================================\n");
        }

        foreach (var s in snps)
        {
            Dist(s, out long dist);
            sb.Append(L(s.PosR, 8)).Append("   ").Append(s.CharR).Append(' ').Append(s.CharQ).Append("   ").Append(L(s.PosQ).PadRight(8)).Append("  | ");
            sb.Append(L(s.Buff, 8)).Append(' ').Append(L(dist, 8)).Append("  | ");
            if (o.ShowConflict)
                sb.Append(L(s.ConflictR, 4)).Append(' ').Append(L(s.ConflictQ, 4)).Append("  | ");
            if (o.ShowLength)
                sb.Append(L(s.Edge.ReferenceNode.Length, 8)).Append(' ').Append(L(s.Edge.QueryNode.Length, 8)).Append("  | ");
            if (o.Context != 0)
            {
                sb.Append(' ', ctxc - ctxw).Append(' ').Append(s.ContextR).Append("  ").Append(' ', ctxc - ctxw).Append(s.ContextQ).Append("  | ");
            }

            sb.Append(L(s.Edgelet.FrameR, 2)).Append(' ').Append(L(s.Edgelet.FrameQ, 2)).Append("  ");
            sb.Append(s.Edge.ReferenceNode.Id).Append('\t').Append(s.Edge.QueryNode.Id).Append('\n');
        }
    }

    private static void Dist(MummerSnp s, out long dist)
    {
        long distR = s.PosR < s.Edge.ReferenceNode.Length - s.PosR + 1 ? s.PosR : s.Edge.ReferenceNode.Length - s.PosR + 1;
        long distQ = s.PosQ < s.Edge.QueryNode.Length - s.PosQ + 1 ? s.PosQ : s.Edge.QueryNode.Length - s.PosQ + 1;
        dist = distR < distQ ? distR : distQ;
    }

    private static void PrintTabular(StringBuilder sb, List<MummerSnp> snps, MummerDeltaGraph graph, MummerShowSnpsOptions o)
    {
        if (o.PrintHeader)
        {
            sb.Append(graph.ReferencePath).Append(' ').Append(graph.QueryPath).Append('\n').Append(graph.DataType).Append("\n\n");
            sb.Append("[P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t");
            if (o.ShowConflict)
                sb.Append("[R]\t[Q]\t");
            if (o.ShowLength)
                sb.Append("[LEN R]\t[LEN Q]\t");
            if (o.Context != 0)
                sb.Append("[CTX R]\t[CTX Q]\t");
            sb.Append("[FRM]\t[TAGS]\n");
        }

        foreach (var s in snps)
        {
            Dist(s, out long dist);
            sb.Append(L(s.PosR)).Append('\t').Append(s.CharR).Append('\t').Append(s.CharQ).Append('\t').Append(L(s.PosQ)).Append('\t');
            sb.Append(L(s.Buff)).Append('\t').Append(L(dist)).Append('\t');
            if (o.ShowConflict)
                sb.Append(L(s.ConflictR)).Append('\t').Append(L(s.ConflictQ)).Append('\t');
            if (o.ShowLength)
                sb.Append(L(s.Edge.ReferenceNode.Length)).Append('\t').Append(L(s.Edge.QueryNode.Length)).Append('\t');
            if (o.Context != 0)
                sb.Append(s.ContextR).Append('\t').Append(s.ContextQ).Append('\t');
            sb.Append(L(s.Edgelet.FrameR)).Append('\t').Append(L(s.Edgelet.FrameQ)).Append('\t');
            sb.Append(s.Edge.ReferenceNode.Id).Append('\t').Append(s.Edge.QueryNode.Id).Append('\n');
        }
    }
}
