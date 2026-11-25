using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>Options of <see cref="MummerShowDiff.Run"/> (<c>show-diff</c> 3.23).</summary>
public sealed record MummerShowDiffOptions
{
    /// <summary><c>-f</c>: AMOS-style output.</summary>
    public bool Amos { get; init; }

    /// <summary><c>-r</c>: the differences of the reference (the default).</summary>
    public bool ReferenceDiff { get; init; }

    /// <summary><c>-q</c>: the differences of the query (cannot be combined with -r).</summary>
    public bool QueryDiff { get; init; }

    /// <summary><c>-H</c> inverted: print the header (default true).</summary>
    public bool PrintHeader { get; init; } = true;
}

/// <summary>
/// <c>show-diff</c> of MUMmer 3.23 (<c>src/tigr/show-diff.cc</c>): the structural differences between the sequences of a delta file - breaks (BRK), gaps (GAP), jumps (JMP), inversions (INV), duplications (DUP)
/// and jumps between sequences (SEQ) - read off the many-to-many mapping (<c>flagMtoM</c>, then <c>clean</c>) of the alignments along every reference (or, with -q, query) record.
/// </summary>
public static class MummerShowDiff
{
    /// <summary>Runs <c>show-diff</c> on the text of a delta file.</summary>
    /// <param name="deltaText">The delta file.</param>
    /// <param name="options">Switches (null = defaults: the reference differences).</param>
    /// <exception cref="ArgumentException">-r and -q together (the program's error).</exception>
    /// <exception cref="FormatException">The delta file cannot be parsed.</exception>
    public static MummerTextOutput Run(string deltaText, MummerShowDiffOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(deltaText);
        var o = options ?? new MummerShowDiffOptions();
        if (o.ReferenceDiff && o.QueryDiff)
            throw new ArgumentException("ERROR: -r and -q options cannot be combined", nameof(options));
        bool refDiff = o.ReferenceDiff || !o.QueryDiff;
        bool qryDiff = o.QueryDiff;

        var graph = MummerDeltaGraph.Build(deltaText, false);
        graph.FlagMtoM();
        graph.Clean();

        var sb = new StringBuilder();
        if (o.PrintHeader)
        {
            sb.Append(graph.ReferencePath).Append(' ').Append(graph.QueryPath).Append('\n').Append(graph.DataType).Append("\n\n");
            sb.Append("[SEQ]\t[TYPE]\t[S1]\t[E1]\t[LEN 1]\n");
        }

        if (refDiff)
            Diff(sb, graph, o.Amos, true);
        if (qryDiff)
            Diff(sb, graph, o.Amos, false);
        return new MummerTextOutput(sb.ToString(), string.Empty);
    }

    private static bool EdgeletLess(MummerDeltaEdgelet i, MummerDeltaEdgelet j, bool byRef)
    {
        // EdgeletIdQLoQCmp_t / EdgeletIdRLoRCmp_t: the nodes are compared by address (allocation order); an alignment without edge (the padding) comes first
        if (i.Edge != null && j.Edge != null)
        {
            int a = byRef ? i.Edge.ReferenceNode.Order : i.Edge.QueryNode.Order;
            int b = byRef ? j.Edge.ReferenceNode.Order : j.Edge.QueryNode.Order;
            if (a < b)
                return true;
            if (a > b)
                return false;
        }
        else if (i.Edge == null && j.Edge != null)
        {
            return true;
        }
        else if (i.Edge != null && j.Edge == null)
        {
            return false;
        }

        return byRef ? i.LoR < j.LoR : i.LoQ < j.LoQ;
    }

    private static void Diff(StringBuilder sb, MummerDeltaGraph graph, bool amos, bool reference)
    {
        var lpad = new MummerDeltaEdgelet { IsRLis = true, IsQLis = true };
        var rpad = new MummerDeltaEdgelet { IsRLis = true, IsQLis = true };
        var nodes = reference ? graph.ReferenceNodes : graph.QueryNodes;
        foreach (var kv in nodes)
        {
            string id = kv.Key;
            var aligns = new List<MummerDeltaEdgelet>();
            foreach (var ei in kv.Value.Edges)
                foreach (var el in ei.Edgelets)
                    aligns.Add(el);

            if (reference)
            {
                rpad.LoR = rpad.HiR = kv.Value.Length + 1;
                rpad.LoQ = rpad.HiQ = long.MaxValue;
            }
            else
            {
                rpad.LoQ = rpad.HiQ = kv.Value.Length + 1;
                rpad.LoR = rpad.HiR = long.MaxValue;
            }

            aligns.Add(lpad);
            aligns.Add(rpad);
            int nAligns = aligns.Count;

            LibstdcxxSort.Sort(aligns, (a, b) => EdgeletLess(a, b, !reference));
            long j = 0;
            foreach (var al in aligns)
            {
                bool lis = reference ? al.IsQLis : al.IsRLis;
                al.StopCount = lis ? j++ : -1;
            }

            LibstdcxxSort.Sort(aligns, reference ? (a, b) => a.LoR < b.LoR : (a, b) => a.LoQ < b.LoQ);
            if (!ReferenceEquals(aligns[0], lpad) || !ReferenceEquals(aligns[nAligns - 1], rpad))
                throw new InvalidOperationException("Assertion `aligns[0] == &lpad && aligns[nAligns-1] == &rpad' failed.");

            var pa = aligns[0];
            var pga = aligns[0];
            for (int i = 1; i != nAligns; ++i)
            {
                var a = aligns[i];
                if (!(reference ? a.IsRLis : a.IsQLis))
                    continue;
                long Lo(MummerDeltaEdgelet x) => reference ? x.LoR : x.LoQ;
                long Hi(MummerDeltaEdgelet x) => reference ? x.HiR : x.HiQ;
                long LoO(MummerDeltaEdgelet x) => reference ? x.LoQ : x.LoR;
                long HiO(MummerDeltaEdgelet x) => reference ? x.HiQ : x.HiR;
                bool otherLis = reference ? a.IsQLis : a.IsRLis;
                long gap = Lo(a) - Hi(pa) - 1;

                if (a.Edge == null)
                {
                    Brk(sb, amos, id, Hi(pa) + 1, Lo(a) - 1);
                }
                else if (otherLis && a.Edge == pga.Edge)
                {
                    if (a.Slope != pga.Slope || a.StopCount != pga.StopCount + pga.Slope)
                    {
                        if (a.Slope == pga.Slope)
                            Simple(sb, amos, id, "JMP", Hi(pa) + 1, Lo(a) - 1);
                        else
                            Simple(sb, amos, id, "INV", Hi(pa) + 1, Lo(a) - 1);
                    }
                    else if (ReferenceEquals(pa, pga))
                    {
                        long gapO = a.IsPositive ? LoO(a) - HiO(pga) - 1 : LoO(pga) - HiO(a) - 1;
                        Gap(sb, amos, id, Hi(pa) + 1, Lo(a) - 1, gap, gapO);
                    }
                    else
                    {
                        Brk(sb, amos, id, Hi(pa) + 1, Lo(a) - 1);
                    }
                }
                else if (!otherLis)
                {
                    Brk(sb, amos, id, Hi(pa) + 1, Lo(a) - 1);
                    Simple(sb, amos, id, "DUP", Lo(a), Hi(a));
                }
                else if (pga.Edge != null)
                {
                    string p = reference ? pga.Edge.QueryNode.Id : pa.Edge!.ReferenceNode.Id;
                    string n = reference ? a.Edge.QueryNode.Id : a.Edge.ReferenceNode.Id;
                    SeqJmp(sb, amos, id, p, n, Hi(pa) + 1, Lo(a) - 1);
                }
                else
                {
                    Brk(sb, amos, id, Hi(pa) + 1, Lo(a) - 1);
                }

                if (otherLis)
                    pga = a;
                pa = a;
            }
        }
    }

    private static string L(long v) => v.ToString(CultureInfo.InvariantCulture);

    private static void Brk(StringBuilder sb, bool amos, string seq, long s, long e)
    {
        if (e - s + 1 <= 0)
            return;
        Simple(sb, amos, seq, "BRK", s, e);
    }

    private static void Simple(StringBuilder sb, bool amos, string seq, string type, long s, long e)
    {
        if (!amos)
            sb.Append(seq).Append('\t').Append(type).Append('\t').Append(L(s)).Append('\t').Append(L(e)).Append('\t').Append(L(e - s + 1)).Append('\n');
        else
            sb.Append(seq).Append("\tA\t").Append(L(s)).Append('\t').Append(L(e)).Append('\t').Append(type).Append('\t').Append(L(s)).Append('\t').Append(L(e)).Append('\t').Append(L(e - s + 1)).Append('\n');
    }

    private static void SeqJmp(StringBuilder sb, bool amos, string seq, string seqp, string seqn, long s, long e)
    {
        if (!amos)
            sb.Append(seq).Append("\tSEQ\t").Append(L(s)).Append('\t').Append(L(e)).Append('\t').Append(L(e - s + 1)).Append('\t').Append(seqp).Append('\t').Append(seqn).Append('\n');
        else
            sb.Append(seq).Append("\tA\t").Append(L(s)).Append('\t').Append(L(e)).Append("\tSEQ\t").Append(L(s)).Append('\t').Append(L(e)).Append('\t').Append(L(e - s + 1)).Append('\t').Append(seqp).Append('\t').Append(seqn).Append('\n');
    }

    public static void Gap(StringBuilder sb, bool amos, string seq, long s, long e, long gap1, long gap2)
    {
        if (!amos)
            sb.Append(seq).Append("\tGAP\t").Append(L(s)).Append('\t').Append(L(e)).Append('\t').Append(L(gap1)).Append('\t').Append(L(gap2)).Append('\t').Append(L(gap1 - gap2)).Append('\n');
        else
            sb.Append(seq).Append("\tA\t").Append(L(s)).Append('\t').Append(L(e)).Append("\tGAP\t").Append(L(s)).Append('\t').Append(L(e)).Append('\t').Append(L(gap1)).Append('\t').Append(L(gap2)).Append('\t').Append(L(gap1 - gap2)).Append('\n');
    }
}
