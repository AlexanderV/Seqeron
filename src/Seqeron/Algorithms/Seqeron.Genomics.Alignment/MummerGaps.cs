using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// <c>gaps</c> of MUMmer 3.23 (<c>src/tigr/gaps.cc</c>, Delcher): reads a list of matches between two strings (<c>start1 start2 length</c>, as <c>mummer</c> prints them without the header lines) and prints first the longest consistent
/// (ordered along both strings) set of matches and then all the other matches, in the format of <c>mgaps</c> (<c>start1 start2 length overlap gap1 gap2</c>).
/// </summary>
public static class MummerGaps
{
    private sealed class Match
    {
        public long Start1, Start2, Len;
        public long SimpleScore, SimpleFrom, SimpleAdj;
        public bool Good;
    }

    /// <summary>Runs <c>gaps &lt;genomefile&gt; [-r] &lt; matches</c>.</summary>
    /// <param name="fileName">The first argument of the program, printed in the header lines.</param>
    /// <param name="reverseComplement">The second argument is <c>-r</c> (the header lines read <c>reverse</c>).</param>
    /// <param name="matches">The matches, white-space separated numbers three by three.</param>
    /// <exception cref="FormatException">The input is not a list of numbers (the program never terminates on it) or empty (the program crashes with a segmentation fault).</exception>
    public static string Run(string fileName, bool reverseComplement, string matches)
    {
        ArgumentNullException.ThrowIfNull(fileName);
        ArgumentNullException.ThrowIfNull(matches);
        var a = new List<Match>();
        var tokens = matches.Split([' ', '\t', '\n', '\r', '\v', '\f'], StringSplitOptions.RemoveEmptyEntries);
        long s1 = 0, s2 = 0, len = 0;
        int p = 0;
        while (p < tokens.Length)
        {
            // scanf ("%ld %ld %ld"): a partial read leaves the earlier values in place and still counts as an item
            int read = 0;
            while (read < 3 && p < tokens.Length)
            {
                if (!long.TryParse(tokens[p], NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out long v))
                    throw new FormatException("gaps: the input contains '" + tokens[p] + "': the program loops forever on text that is not a number.");
                p++;
                if (read == 0)
                    s1 = v;
                else if (read == 1)
                    s2 = v;
                else
                    len = v;
                read++;
            }

            a.Add(new Match { Start1 = s1, Start2 = s2, Len = len });
        }

        if (a.Count == 0)
            throw new FormatException("gaps: no matches (the program dereferences a null array and crashes).");

        // qsort by Start2: glibc's qsort is a merge sort (stable)
        var m = a.Select((x, i) => (x, i)).OrderBy(t => t.x.Start2).ThenBy(t => t.i).Select(t => t.x).ToArray();
        for (int i = 0; i < m.Length; i++)
        {
            m[i].SimpleScore = m[i].Len;
            m[i].SimpleAdj = 0;
            m[i].SimpleFrom = -1;
            for (int j = 0; j < i; j++)
            {
                long olap2 = m[j].Start2 + m[j].Len - m[i].Start2;
                long olap = Math.Max(0, olap2);
                long olap1 = m[j].Start1 + m[j].Len - m[i].Start1;
                olap = Math.Max(olap, olap1);
                if (m[j].SimpleScore + m[i].Len - olap > m[i].SimpleScore)
                {
                    m[i].SimpleFrom = j;
                    m[i].SimpleScore = m[j].SimpleScore + m[i].Len - olap;
                    m[i].SimpleAdj = olap;
                }
            }
        }

        int best = 0;
        for (int i = 1; i < m.Length; i++)
            if (m[i].SimpleScore > m[best].SimpleScore)
                best = i;
        for (long i = best; i >= 0; i = m[i].SimpleFrom)
            m[i].Good = true;

        var ci = CultureInfo.InvariantCulture;
        var sb = new StringBuilder();
        string tag = reverseComplement ? "reverse" : string.Empty;
        sb.Append("> ").Append(fileName).Append(' ').Append(tag).Append(" Consistent matches\n");
        long prev = -1;
        for (int i = 0; i < m.Length; i++)
        {
            if (!m[i].Good)
                continue;
            if (prev == -1)
                sb.Append(m[i].Start1.ToString(ci).PadLeft(8)).Append(' ').Append(m[i].Start2.ToString(ci).PadLeft(8)).Append(' ').Append(m[i].Len.ToString(ci).PadLeft(6)).Append(' ').Append("none".PadLeft(7)).Append(' ').Append("-".PadLeft(6)).Append(' ').Append("-".PadLeft(6)).Append('\n');
            else
            {
                long adj = m[i].SimpleAdj;
                AppendJoined(sb, m[i], adj, m[prev]);
            }

            prev = i;
        }

        sb.Append("> ").Append(fileName).Append(' ').Append(tag).Append(" Other matches\n");
        prev = -1;
        for (int i = 0; i < m.Length; i++)
        {
            if (m[i].Good)
                continue;
            if (prev == -1)
                sb.Append(m[i].Start1.ToString(ci).PadLeft(8)).Append(' ').Append(m[i].Start2.ToString(ci).PadLeft(8)).Append(' ').Append(m[i].Len.ToString(ci).PadLeft(6)).Append(' ').Append("none".PadLeft(7)).Append(' ').Append("-".PadLeft(6)).Append(' ').Append("-".PadLeft(6)).Append('\n');
            else
            {
                bool joined = m[i].SimpleFrom == prev;
                long adj = joined ? m[i].SimpleAdj : 0;
                if (joined)
                    AppendJoined(sb, m[i], adj, m[prev]);
                else
                {
                    sb.Append((m[i].Start1 + adj).ToString(ci).PadLeft(8)).Append(' ').Append((m[i].Start2 + adj).ToString(ci).PadLeft(8)).Append(' ').Append((m[i].Len - adj).ToString(ci).PadLeft(6))
                        .Append(' ').Append("none".PadLeft(7)).Append(' ').Append("-".PadLeft(6)).Append(' ').Append("-".PadLeft(6)).Append('\n');
                }
            }

            prev = i;
        }

        return sb.ToString();

        static void AppendJoined(StringBuilder sb, Match cur, long adj, Match before)
        {
            var ci = CultureInfo.InvariantCulture;
            sb.Append((cur.Start1 + adj).ToString(ci).PadLeft(8)).Append(' ').Append((cur.Start2 + adj).ToString(ci).PadLeft(8)).Append(' ').Append((cur.Len - adj).ToString(ci).PadLeft(6))
                .Append(' ').Append((adj == 0 ? "none" : (-adj).ToString(ci)).PadLeft(7)).Append(' ')
                .Append((cur.Start1 + adj - before.Start1 - before.Len).ToString(ci).PadLeft(6)).Append(' ')
                .Append((cur.Start2 + adj - before.Start2 - before.Len).ToString(ci).PadLeft(6)).Append('\n');
        }
    }
}

/// <summary>
/// <c>delta2blocks</c> of the Debian MUMmer 3.23 package (<c>addition_from_mugsy.patch</c>, <c>src/tigr/delta2blocks.cc</c>): the alignments of a nucmer delta file between one reference and one query record as one line each,
/// <c>refId qryId sR eR sQ eQ CIGAR</c> with the operations M (matches), X (mismatches), I (an extra character in the reference, a positive delta) and D (an extra character in the query, a negative delta).
/// </summary>
public static class MummerDelta2Blocks
{
    /// <summary>Runs <c>delta2blocks [-q|-r] [-w int] [-x int] deltafile refId qryId</c>.</summary>
    /// <exception cref="FormatException">The delta file cannot be parsed or an id is missing.</exception>
    /// <exception cref="InvalidOperationException">The delta file is a promer file: the program dereferences a null pointer (SIGSEGV).</exception>
    public static MummerTextOutput Run(string deltaText, string referenceFasta, string queryFasta, string referenceId, string queryId, MummerShowAlignsOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(deltaText);
        ArgumentNullException.ThrowIfNull(referenceFasta);
        ArgumentNullException.ThrowIfNull(queryFasta);
        var o = options ?? new MummerShowAlignsOptions();
        var err = new StringBuilder();
        if (o.ScreenWidth <= 11)
            err.Append("WARNING: invalid screen width ").Append(o.ScreenWidth.ToString(CultureInfo.InvariantCulture)).Append(", using default\n");
        if (o.MatrixType < 1 || o.MatrixType > 3)
            err.Append("WARNING: invalid matrix type ").Append(o.MatrixType.ToString(CultureInfo.InvariantCulture)).Append(", using default\n");
        if (o.SortByQuery && o.SortByReference)
            err.Append("WARNING: both -r and -q were passed, -q ignored\n");

        var file = MummerDeltaFile.Parse(deltaText, true);
        var found = file.Records.FirstOrDefault(r => r.ReferenceId == referenceId && r.QueryId == queryId)
            ?? throw new FormatException("ERROR: Could not find any alignments for " + referenceId + " and " + queryId);
        var aligns = found.Alignments.Select(x => (SR: x.StartR, ER: x.EndR, SQ: x.StartQ, EQ: x.EndQ, Delta: x.Deltas)).ToList();
        byte[] r = Find(referenceFasta, referenceId, "reference");
        byte[] q = Find(queryFasta, queryId, "query");
        if (file.DataType != "NUCMER")
            throw new InvalidOperationException("delta2blocks: the program reads the sequences of a promer alignment through a null pointer and crashes (SIGSEGV).");
        if (o.SortByReference)
            LibstdcxxSort.Sort(aligns, (a, b) => a.SR < b.SR);
        else if (o.SortByQuery)
            LibstdcxxSort.Sort(aligns, (a, b) => a.SQ < b.SQ);

        var sb = new StringBuilder();
        var ci = CultureInfo.InvariantCulture;
        long lenR = r.Length - 2, lenQ = q.Length - 2;
        var a1 = r;
        var a4 = MummerSequence.ReverseComplement(r);
        var b1 = q;
        var b4 = MummerSequence.ReverseComplement(q);
        foreach (var ap in aligns)
        {
            // printf ("%s %s %d %d %d %d ", ...): %d reads the low 32 bits of the long values
            sb.Append(referenceId).Append(' ').Append(queryId).Append(' ').Append(((int)ap.SR).ToString(ci)).Append(' ').Append(((int)ap.ER).ToString(ci)).Append(' ')
                .Append(((int)ap.SQ).ToString(ci)).Append(' ').Append(((int)ap.EQ).ToString(ci)).Append(' ');
            long sR = ap.SR, eR = ap.ER, sQ = ap.SQ;
            bool revR = sR > eR, revQ = ap.SQ > ap.EQ;
            if (revR)
            {
                sR = lenR - sR + 1;
                eR = lenR - eR + 1;
            }

            if (revQ)
                sQ = lenQ - sQ + 1;
            byte[] aa = revR ? a4 : a1, bb = revQ ? b4 : b1;
            long apos = sR, bpos = sQ, remain = eR - sR + 1;
            char kind = '\0';
            long count = 0;
            void Op(char k)
            {
                if (kind == k)
                {
                    count++;
                    return;
                }

                if (count > 0)
                    sb.Append(count.ToString(ci)).Append(kind);
                kind = k;
                count = 1;
            }

            foreach (long dp in ap.Delta)
            {
                if (dp == 0)
                    break;
                long delta = Math.Abs(dp);
                long i;
                for (i = 1; i < delta; i++)
                {
                    Op(aa[apos] == bb[bpos] ? 'M' : 'X');
                    apos++;
                    bpos++;
                }

                remain -= i - 1;
                if (dp > 0)
                {
                    Op('I');
                    apos++;
                    remain--;
                }
                else
                {
                    Op('D');
                    bpos++;
                }
            }

            for (long i = 0; i < remain; i++)
            {
                Op(aa[apos] == bb[bpos] ? 'M' : 'X');
                apos++;
                bpos++;
            }

            if (count > 0)
                sb.Append(count.ToString(ci)).Append(kind);
            sb.Append('\n');
        }

        return new MummerTextOutput(sb.ToString(), err.ToString());
    }

    private static byte[] Find(string fasta, string id, string what)
    {
        foreach (var (rid, seq) in NucmerPostnuc.ReadFasta(fasta, out _))
        {
            if (rid == id)
            {
                var bytes = new byte[seq.Length + 2];
                for (int i = 0; i < seq.Length; i++)
                    bytes[i + 1] = (byte)seq[i];
                return bytes;
            }
        }

        throw new FormatException("ERROR: Could not find " + id + " in the " + what + " file");
    }
}
