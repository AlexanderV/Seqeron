using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>The switches of <see cref="MummerDelta2Maf.Run"/> (<c>delta2maf</c> of the Debian MUMmer 3.23 package).</summary>
public sealed record MummerDelta2MafOptions
{
    /// <summary><c>-q</c>: sort the alignments by the query start (ignored with <c>-r</c>).</summary>
    public bool SortByQuery { get; init; }

    /// <summary><c>-r</c>: sort the alignments by the reference start.</summary>
    public bool SortByReference { get; init; }

    /// <summary><c>-w</c>: the buffer width (default 100000; a value of 11 or less is replaced with a warning); it does not change the output.</summary>
    public int ScreenWidth { get; init; } = 100000;

    /// <summary><c>-x</c>: the matrix type (1 - 3, only checked: the alignments are always written as DNA).</summary>
    public int MatrixType { get; init; } = 2;
}

/// <summary>
/// <c>delta2maf</c> of the Debian MUMmer 3.23 package (<c>addition_from_mugsy.patch</c>, <c>0009-improve-performance-of-delta2maf.patch</c>, <c>src/tigr/delta2maf.cc</c>): the alignments of a nucmer delta file as a MAF
/// (multiple alignment format) file: <c>##maf version=1 scoring=single_cov2</c>, per alignment <c>a score=n</c> (n = |eR - sR|) and two <c>s</c> lines (source, 0-based start, length, strand, source length, the gapped text), and
/// <c>##eof maf</c>. The reference and query files named in the delta header are read into one map by id (the first record of an id wins: a query record with the id of a reference record is read as that reference record).
/// </summary>
public static class MummerDelta2Maf
{
    private sealed class Aln
    {
        public long SQ, EQ, SR, ER;
        public List<long> Delta = [];
        public string IdR = string.Empty, IdQ = string.Empty;
    }

    private static byte[] Pick(bool reverse, byte[] reversed, byte[] forward) => reverse ? reversed : forward;

    private static byte[] ReverseOf(Dictionary<string, byte[]> cache, string id, byte[] sequence)
    {
        if (!cache.TryGetValue(id, out byte[]? rc))
        {
            rc = MummerSequence.ReverseComplement(sequence);
            cache[id] = rc;
        }

        return rc;
    }

    /// <summary>Runs <c>delta2maf [options] deltafile [refId] [qryId]</c>.</summary>
    /// <param name="deltaText">The delta file (nucmer).</param>
    /// <param name="referenceFasta">The reference FASTA named in the delta header.</param>
    /// <param name="queryFasta">The query FASTA named in the delta header.</param>
    /// <param name="referenceId">Only the alignments of this reference record (empty = any).</param>
    /// <param name="queryId">Only the alignments of this query record (empty = any; with a reference id both must match).</param>
    /// <param name="options">The switches.</param>
    /// <exception cref="FormatException">The delta file cannot be parsed, no alignment is found, or a record named in the delta file is not in the FASTA files (assertion).</exception>
    /// <exception cref="NotSupportedException">The delta file is a promer file: the program reads its amino-acid deltas as nucleotide deltas and runs past the sequences (undefined).</exception>
    public static MummerTextOutput Run(string deltaText, string referenceFasta, string queryFasta, string referenceId = "", string queryId = "", MummerDelta2MafOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(deltaText);
        ArgumentNullException.ThrowIfNull(referenceFasta);
        ArgumentNullException.ThrowIfNull(queryFasta);
        ArgumentNullException.ThrowIfNull(referenceId);
        ArgumentNullException.ThrowIfNull(queryId);
        var o = options ?? new MummerDelta2MafOptions();
        var err = new StringBuilder();
        var ci = CultureInfo.InvariantCulture;
        if (o.ScreenWidth <= 11)
            err.Append("WARNING: invalid screen width ").Append("100000").Append(", using default\n");
        if (o.MatrixType < 1 || o.MatrixType > 3)
            err.Append("WARNING: invalid matrix type ").Append(o.MatrixType.ToString(ci)).Append(", using default\n");
        if (o.SortByQuery && o.SortByReference)
            err.Append("WARNING: both -r and -q were passed, -q ignored\n");

        var file = MummerDeltaFile.Parse(deltaText, true);
        if (file.DataType != "NUCMER")
            throw new NotSupportedException("delta2maf of a promer delta file: the program reads the amino-acid deltas as nucleotide deltas and runs past the sequences.");
        var aligns = new List<Aln>();
        foreach (var rec in file.Records)
        {
            bool found;
            if (referenceId.Length == 0 && queryId.Length == 0)
                found = true;
            else if (referenceId.Length != 0 && queryId.Length != 0)
                found = rec.ReferenceId == referenceId && rec.QueryId == queryId;
            else
                found = referenceId.Length != 0 ? rec.ReferenceId == referenceId : rec.QueryId == queryId;
            if (!found)
                continue;
            foreach (var a in rec.Alignments)
                aligns.Add(new Aln { SR = a.StartR, ER = a.EndR, SQ = a.StartQ, EQ = a.EndQ, Delta = a.Deltas, IdR = rec.ReferenceId, IdQ = rec.QueryId });
        }

        if (aligns.Count == 0)
            throw new FormatException("ERROR: Could not find any alignments for " + referenceId + " and " + queryId);

        var seqs = new Dictionary<string, byte[]>(StringComparer.Ordinal);
        foreach (string fasta in new[] { referenceFasta, queryFasta })
        {
            foreach (var (id, seq) in NucmerPostnuc.ReadFasta(fasta, out _))
            {
                if (seqs.ContainsKey(id))
                    continue;   // map::insert keeps the first record of an id
                var bytes = new byte[seq.Length + 2];
                for (int i = 0; i < seq.Length; i++)
                    bytes[i + 1] = (byte)seq[i];
                seqs[id] = bytes;
            }
        }

        if (o.SortByReference)
            LibstdcxxSort.Sort(aligns, (a, b) => a.SR < b.SR);
        else if (o.SortByQuery)
            LibstdcxxSort.Sort(aligns, (a, b) => a.SQ < b.SQ);

        var sb = new StringBuilder("##maf version=1 scoring=single_cov2\n");
        var reverse = new Dictionary<string, byte[]>(StringComparer.Ordinal);
        foreach (var ap in aligns)
        {
            if (!seqs.TryGetValue(ap.IdR, out byte[]? r) || !seqs.TryGetValue(ap.IdQ, out byte[]? q))
                throw new FormatException("Assertion `finditer != seqsMap.end()' failed.");
            long lenR = r.Length - 2, lenQ = q.Length - 2;
            byte[] r4 = ReverseOf(reverse, ap.IdR, r);
            byte[] q4 = ReverseOf(reverse, ap.IdQ, q);
            long sR = ap.SR, eR = ap.ER, sQ = ap.SQ, eQ = ap.EQ;
            bool negR = false, negQ = false;
            if (sR > eR)
            {
                sR = lenR - sR + 1;
                eR = lenR - eR + 1;
                negR = true;
            }

            if (sQ > eQ)
            {
                sQ = lenQ - sQ + 1;
                eQ = lenQ - eQ + 1;
                negQ = true;
            }

            sb.Append("a score=").Append(((int)Math.Abs(ap.ER - ap.SR)).ToString(ci)).Append('\n');
            for (int query = 0; query < 2; query++)
            {
                byte[] seq = query == 0 ? Pick(negR, r4, r) : Pick(negQ, q4, q);
                long pos = query == 0 ? sR : sQ;
                long remain = eR - sR + 1;
                if (query == 1)
                    sb.Append("s ").Append(ap.IdQ).Append(' ').Append(((int)(sQ - 1)).ToString(ci)).Append(' ').Append(((int)(eQ - sQ + 1)).ToString(ci)).Append(' ').Append(negQ ? '-' : '+').Append(' ').Append(((int)lenQ).ToString(ci)).Append(' ');
                else
                    sb.Append("s ").Append(ap.IdR).Append(' ').Append(((int)(sR - 1)).ToString(ci)).Append(' ').Append(((int)(eR - sR + 1)).ToString(ci)).Append(' ').Append(negR ? '-' : '+').Append(' ').Append(((int)lenR).ToString(ci)).Append(' ');
                foreach (long dp in ap.Delta)
                {
                    if (dp == 0)
                        break;
                    long delta = Math.Abs(dp);
                    for (long i = 1; i < delta; i++)
                        sb.Append((char)seq[pos++]);
                    remain -= delta - 1;
                    if (dp > 0)
                    {
                        sb.Append(query == 0 ? (char)seq[pos++] : '-');
                        remain--;
                    }
                    else
                        sb.Append(query == 0 ? '-' : (char)seq[pos++]);
                }

                for (long i = 0; i < remain; i++)
                    sb.Append((char)seq[pos++]);
                sb.Append('\n');
                if (query == 1)
                    sb.Append('\n');
            }
        }

        sb.Append("##eof maf");
        return new MummerTextOutput(sb.ToString(), err.ToString());
    }
}
