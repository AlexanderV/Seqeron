using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>Options of <see cref="MummerShowAligns.Run"/> (<c>show-aligns</c> 3.23).</summary>
public sealed record MummerShowAlignsOptions
{
    /// <summary><c>-q</c>: sort the alignments by the query start.</summary>
    public bool SortByQuery { get; init; }

    /// <summary><c>-r</c>: sort the alignments by the reference start.</summary>
    public bool SortByReference { get; init; }

    /// <summary><c>-w</c>: screen width (default 60; a value of 11 or less is replaced by 60 with a warning).</summary>
    public int ScreenWidth { get; init; } = 60;
}

/// <summary>
/// <c>show-aligns</c> of MUMmer 3.23 (<c>src/tigr/show-aligns.cc</c>): the alignments between one reference record and one query record of a nucmer delta file, drawn as three rows (reference, query, a
/// <c>^</c> under every difference; an indel is a <c>.</c> in the other row) in blocks of the screen width, each block row prefixed with the 1-based position of its first base on the forward strand.
/// </summary>
public static class MummerShowAligns
{
    private const int LinePrefixLen = 11;

    private sealed class Aln
    {
        public long SQ, EQ, SR, ER;
        public List<long> Delta = [];
    }

    /// <summary>Runs <c>show-aligns</c>.</summary>
    /// <param name="deltaText">The delta file (nucmer).</param>
    /// <param name="referenceFasta">The reference FASTA named in the delta header.</param>
    /// <param name="queryFasta">The query FASTA named in the delta header.</param>
    /// <param name="referenceId">Reference record id.</param>
    /// <param name="queryId">Query record id.</param>
    /// <param name="options">Switches (null = defaults).</param>
    /// <exception cref="FormatException">The delta file cannot be parsed or an id is missing (the program's error exits).</exception>
    /// <exception cref="NotSupportedException">The delta file is a promer file (the amino-acid alignment is a registered sub-item).</exception>
    public static MummerTextOutput Run(string deltaText, string referenceFasta, string queryFasta, string referenceId, string queryId, MummerShowAlignsOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(deltaText);
        ArgumentNullException.ThrowIfNull(referenceFasta);
        ArgumentNullException.ThrowIfNull(queryFasta);
        var o = options ?? new MummerShowAlignsOptions();
        var err = new StringBuilder();
        int width = o.ScreenWidth;
        if (width <= LinePrefixLen)
        {
            err.Append("WARNING: invalid screen width ").Append(width.ToString(CultureInfo.InvariantCulture)).Append(", using default\n");
            width = 60;
        }

        if (o.SortByQuery && o.SortByReference)
            err.Append("WARNING: both -r and -q were passed, -q ignored\n");

        var file = MummerDeltaFile.Parse(deltaText, true);
        if (file.DataType != "NUCMER")
            throw new NotSupportedException("show-aligns of a promer delta file (amino-acid alignment) is not ported.");
        MummerDeltaRecord? found = file.Records.FirstOrDefault(r => r.ReferenceId == referenceId && r.QueryId == queryId)
            ?? throw new FormatException("ERROR: Could not find any alignments for " + referenceId + " and " + queryId);
        var aligns = found.Alignments.Select(a => new Aln { SR = a.StartR, ER = a.EndR, SQ = a.StartQ, EQ = a.EndQ, Delta = a.Deltas }).ToList();

        byte[] R = Find(referenceFasta, referenceId, "reference");
        byte[] Q = Find(queryFasta, queryId, "query");

        if (o.SortByReference)
            LibstdcxxSort.Sort(aligns, (a, b) => a.SR < b.SR);
        else if (o.SortByQuery)
            LibstdcxxSort.Sort(aligns, (a, b) => a.SQ < b.SQ);

        var sb = new StringBuilder();
        sb.Append(file.ReferencePath).Append(' ').Append(file.QueryPath).Append("\n\n").Append('=', width).Append("\n-- Alignments between ").Append(referenceId).Append(" and ").Append(queryId).Append("\n\n");
        PrintAlignments(sb, aligns, R, Q, width);
        sb.Append('\n').Append('=', width).Append('\n');
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

    private static byte Complement(byte ch)
    {
        switch (char.ToLowerInvariant((char)ch))
        {
            case 'a': return (byte)'t';
            case 'c': return (byte)'g';
            case 'g': return (byte)'c';
            case 't': return (byte)'a';
            case 'r': return (byte)'y';
            case 'y': return (byte)'r';
            case 's': return (byte)'s';
            case 'w': return (byte)'w';
            case 'm': return (byte)'k';
            case 'k': return (byte)'m';
            case 'b': return (byte)'v';
            case 'd': return (byte)'h';
            case 'h': return (byte)'d';
            case 'v': return (byte)'b';
            default: return (byte)char.ToLowerInvariant((char)ch);
        }
    }

    private static byte[] ReverseComplement(byte[] s)
    {
        long len = s.Length - 2;
        var rev = (byte[])s.Clone();
        long lo = 1, hi = len;
        while (lo <= hi)
        {
            byte ch = rev[hi];
            rev[hi] = Complement(rev[lo]);
            rev[lo] = Complement(ch);
            lo++;
            hi--;
        }

        return rev;
    }

    private static long RevC(long coord, long len) => len - coord + 1;

    private static long ToFwd(long coord, long len, int frame) => frame < 0 ? RevC(coord, len) : coord;

    private static void PrintAlignments(StringBuilder sb, List<Aln> aligns, byte[] r, byte[] q, int width)
    {
        long seqLenR = r.Length - 2, seqLenQ = q.Length - 2;
        var a = new byte[7][];
        var b = new byte[7][];
        a[1] = r;
        a[4] = ReverseComplement(r);
        b[1] = q;
        b[4] = ReverseComplement(q);

        var buff1 = new char[width + 2];
        var buff2 = new char[width + 2];
        var buff3 = new char[width + 2];
        for (int i = 0; i < LinePrefixLen; i++)
            buff3[i] = ' ';

        foreach (var ap in aligns)
        {
            long sR = ap.SR, eR = ap.ER, sQ = ap.SQ;
            int frameR = 1, frameQ = 1;
            if (sR > eR)
            {
                sR = RevC(sR, seqLenR);
                eR = RevC(eR, seqLenR);
                frameR += 3;
            }

            if (sQ > ap.EQ)
            {
                sQ = RevC(sQ, seqLenQ);
                frameQ += 3;
            }

            int ai = frameR, bi = frameQ;
            if (frameR > 3)
                frameR = -(frameR - 3);
            if (frameQ > 3)
                frameQ = -(frameQ - 3);

            string tag = string.Concat(frameR > 0 ? "+" : "-", Math.Abs(frameR).ToString(CultureInfo.InvariantCulture), " ", ap.SR.ToString(CultureInfo.InvariantCulture), " - ", ap.ER.ToString(CultureInfo.InvariantCulture),
                " | ", frameQ > 0 ? "+" : "-", Math.Abs(frameQ).ToString(CultureInfo.InvariantCulture), " ", ap.SQ.ToString(CultureInfo.InvariantCulture), " - ", ap.EQ.ToString(CultureInfo.InvariantCulture));
            sb.Append("-- BEGIN alignment [ ").Append(tag).Append(" ]\n\n\n");

            long apos = sR, bpos = sQ;
            long remain = eR - sR + 1;
            int pos;

            void Prefix()
            {
                Put(buff1, (ToFwd(apos, seqLenR, frameR)).ToString(CultureInfo.InvariantCulture).PadRight(10) + " ");
                Put(buff2, (ToFwd(bpos, seqLenQ, frameQ)).ToString(CultureInfo.InvariantCulture).PadRight(10) + " ");
            }

            void Flush(int p)
            {
                buff1[p] = buff2[p] = buff3[p] = '\0';
                sb.Append(Str(buff1)).Append('\n').Append(Str(buff2)).Append('\n').Append(Str(buff3)).Append("\n\n");
                Prefix();
            }

            Prefix();
            pos = LinePrefixLen;

            foreach (long dp in ap.Delta)
            {
                if (dp == 0)
                    break;
                long delta = dp;
                int sign = delta > 0 ? 1 : -1;
                delta = Math.Abs(delta);
                long i;
                for (i = 1; i < delta; i++)
                {
                    if (pos >= width)
                    {
                        Flush(pos);
                        pos = LinePrefixLen;
                    }

                    buff3[pos] = a[ai][apos] == b[bi][bpos] ? ' ' : '^';
                    buff1[pos] = (char)a[ai][apos++];
                    buff2[pos++] = (char)b[bi][bpos++];
                }

                remain -= i - 1;
                if (pos >= width)
                {
                    Flush(pos);
                    pos = LinePrefixLen;
                }

                if (sign == 1)
                {
                    buff3[pos] = '^';
                    buff1[pos] = (char)a[ai][apos++];
                    buff2[pos++] = '.';
                    remain--;
                }
                else
                {
                    buff3[pos] = '^';
                    buff1[pos] = '.';
                    buff2[pos++] = (char)b[bi][bpos++];
                }
            }

            for (long i = 0; i < remain; i++)
            {
                if (pos >= width)
                {
                    Flush(pos);
                    pos = LinePrefixLen;
                }

                buff3[pos] = a[ai][apos] == b[bi][bpos] ? ' ' : '^';
                buff1[pos] = (char)a[ai][apos++];
                buff2[pos++] = (char)b[bi][bpos++];
            }

            if (pos > LinePrefixLen)
                Flush(pos);

            sb.Append("\n--   END alignment [ ").Append(tag).Append(" ]\n");
        }

        static void Put(char[] dst, string s)
        {
            for (int i = 0; i < s.Length && i < dst.Length; i++)
                dst[i] = s[i];
        }

        static string Str(char[] c)
        {
            int n = Array.IndexOf(c, '\0');
            return new string(c, 0, n < 0 ? c.Length : n);
        }
    }
}
