using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>Options of <see cref="MummerShowTiling.Run"/> (<c>show-tiling</c> 3.23). A null value takes the default of the delta file's type (nucmer / promer).</summary>
public sealed record MummerShowTilingOptions
{
    /// <summary><c>-a</c>: print the tab-delimited alignment region coordinates of the tiling path.</summary>
    public bool PrintAlignments { get; init; }

    /// <summary><c>-c</c>: circular reference sequences (tiled contigs may span the origin).</summary>
    public bool Circular { get; init; }

    /// <summary><c>-g</c>: maximum gap between clustered alignments, -1 = infinity (nucmer default 1000, promer -1).</summary>
    public long? MaxGap { get; init; }

    /// <summary><c>-i</c>: minimum percent identity to tile, 0 - 100 (nucmer default 90, promer 55).</summary>
    public float? MinIdentity { get; init; }

    /// <summary><c>-l</c>: minimum contig length to report, -1 = infinity (default 1).</summary>
    public long MinContigLength { get; init; } = 1;

    /// <summary><c>-R</c>: place repetitive contigs at a random one of their copies (implies -V 0).</summary>
    public bool RandomRepeats { get; init; }

    /// <summary><c>-v</c>: minimum contig coverage to tile, 0 - 100 (nucmer default 95, promer 50).</summary>
    public float? MinCoverage { get; init; }

    /// <summary><c>-V</c>: minimum contig coverage difference, 0 - 100 (nucmer default 10, promer 30).</summary>
    public float? MinCoverageDifference { get; init; }

    /// <summary><c>-x</c>: print the XML contig linking information instead of the path.</summary>
    public bool Xml { get; init; }

    /// <summary><c>-t</c>: produce the TIGR-style contig list (the text is returned in <see cref="MummerShowTilingResult.ContigList"/>).</summary>
    public bool ContigList { get; init; }

    /// <summary><c>-u</c>: produce the alignment coordinates of the unusable contigs (<see cref="MummerShowTilingResult.Unusable"/>).</summary>
    public bool Unusable { get; init; }

    /// <summary><c>-p</c>: produce the pseudo molecule of the query contigs (needs the query FASTA; <see cref="MummerShowTilingResult.PseudoMolecule"/>).</summary>
    public bool PseudoMolecule { get; init; }

    /// <summary>The <c>srand (time (NULL))</c> seed of <c>-R</c>; null = the current time.</summary>
    public uint? RandomSeed { get; init; }

    /// <summary>The time of the XML header (<c>ctime</c>); null = now.</summary>
    public DateTime? Now { get; init; }

    /// <summary>The command line words printed in the XML <c>PARAMETERS</c> attribute (argv).</summary>
    public IReadOnlyList<string>? CommandLine { get; init; }
}

/// <summary>Output of <see cref="MummerShowTiling.Run"/>.</summary>
/// <param name="Text">Standard output.</param>
/// <param name="Messages">Standard error.</param>
/// <param name="ContigList">The <c>-t</c> file (null when not requested).</param>
/// <param name="Unusable">The <c>-u</c> file.</param>
/// <param name="PseudoMolecule">The <c>-p</c> file.</param>
public sealed record MummerShowTilingResult(string Text, string Messages, string? ContigList, string? Unusable, string? PseudoMolecule);

/// <summary>
/// <c>show-tiling</c> of MUMmer 3.23 (<c>src/tigr/show-tiling.cc</c>): places the aligned query contigs on the reference by the longest consistent subset of their alignments per reference
/// (a contig is tiled when its coverage is at least <c>-v</c>, better than the second best by <c>-V</c> and its identity at least <c>-i</c>), links the tiled contigs along every reference and prints the tiling path,
/// the alignments (-a), the XML links (-x), a TIGR contig list (-t), the unusable contigs (-u) or the pseudo molecule (-p).
/// </summary>
public static class MummerShowTiling
{
    private const char ForwardChar = '+';
    private const char ReverseChar = '-';
    private const float RepeatPidyDiff = 0.25f;
    private const int UnusableTileLevel = -2;
    private const int IgnoreTileLevel = -1;
    private const int UnusedTileLevel = 0;
    private const int UsedTileLevel = 1;

    private sealed class AlignStats
    {
        public float Idy;
        public char DirQ;
        public long SQ, EQ, SR, ER;
        public string IdR = string.Empty;
        public long SeqLenR;
        public bool IsTiled;
    }

    private sealed class QueryContig
    {
        public string IdQ = string.Empty;
        public long SeqLenQ;
        public byte[]? SeqQ;
        public List<AlignStats> Aligns = [];
        public int TileLevel;
        public char DirQ = '*';
        public string IdR = string.Empty;
        public long SeqLenR, StartR, EndR, LoTrim, HiTrim;
        public QueryContig? LinksTo;
        public bool IsLinkHead;
    }

    private sealed class Ctx
    {
        public bool Circular, RandomRepeats, Promer;
        public long MaxGap, MinContigLength;
        public float MinCoverage, MinCoverageDiff, MinPidy;
        public GlibcRandom Rand = null!;
    }

    private static long RevC(long coord, long len)
    {
        if (len - coord + 1 <= 0)
            throw new InvalidOperationException("Assertion `Len - Coord + 1 > 0' failed.");
        return len - coord + 1;
    }

    private static float Pct(long n, long d)
    {
        float q = (float)n / (float)d;
        return (float)((double)q * 100.0);
    }

    private static float AlignmentQryCoverage(AlignStats a, long seqLenQ) => Pct(a.EQ - a.SQ + 1, seqLenQ);

    private static float SubsetQryCoverage(List<AlignStats> list, int begin, int end, long seqLenQ)
    {
        long cov = 0;
        AlignStats? pre = null;
        for (int i = begin; i < end; i++)
        {
            var ap = list[i];
            if (!ap.IsTiled)
                continue;
            if (pre == null)
            {
                cov += ap.EQ - ap.SQ + 1;
            }
            else
            {
                if (ap.SQ < pre.SQ)
                    throw new InvalidOperationException("Assertion `Ap->sQ >= preAp->sQ' failed.");
                long olap = pre.EQ - ap.SQ + 1;
                olap = olap > 0 ? olap : 0;
                cov += ap.EQ - ap.SQ + 1 - olap;
            }

            pre = ap;
        }

        return Pct(cov, seqLenQ);
    }

    private static float SubsetQrySyntenyCoverage(List<AlignStats> list, int begin, int end, long seqLenQ)
    {
        long low = long.MaxValue, high = -long.MaxValue;
        for (int i = begin; i < end; i++)
        {
            if (!list[i].IsTiled)
                continue;
            if (list[i].SQ < low)
                low = list[i].SQ;
            if (list[i].EQ > high)
                high = list[i].EQ;
        }

        return Pct(high - low + 1, seqLenQ);
    }

    private static float SubsetIdentity(List<AlignStats> list, int begin, int end)
    {
        long n = 0;
        float tot = 0;
        for (int i = begin; i < end; i++)
        {
            var ap = list[i];
            if (!ap.IsTiled)
                continue;
            long len = ap.EQ - ap.SQ + 1;
            n += len;
            float term = ap.Idy * len;
            tot += term;
        }

        return tot / (float)n;
    }

    /// <summary>Runs <c>show-tiling</c>.</summary>
    /// <param name="deltaText">The delta file.</param>
    /// <param name="queryFasta">The query FASTA named in the header (only for <see cref="MummerShowTilingOptions.PseudoMolecule"/>).</param>
    /// <param name="options">Switches (null = defaults).</param>
    /// <exception cref="ArgumentOutOfRangeException">A range error of the program's option check.</exception>
    /// <exception cref="FormatException">The delta file cannot be parsed or a query sequence is missing for -p.</exception>
    public static MummerShowTilingResult Run(string deltaText, string? queryFasta = null, MummerShowTilingOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(deltaText);
        var o = options ?? new MummerShowTilingOptions();
        float minPidy = o.MinIdentity ?? 0, minCov = o.MinCoverage ?? 0, minCovDiff = o.MinCoverageDifference ?? 0;
        if (minPidy < 0.0 || minPidy > 100.0 || minCov < 0.0 || minCov > 100.0 || minCovDiff < 0.0 || minCovDiff > 100.0)
            throw new ArgumentOutOfRangeException(nameof(options), "ERROR: Percentages must be within the range [0.0, 100.0]");
        if (o.MinContigLength < -1 || (o.MaxGap ?? 0) < -1)
            throw new ArgumentOutOfRangeException(nameof(options), "ERROR: Size values must be within the range [-1, oo]");
        var cx = new Ctx { Circular = o.Circular, RandomRepeats = o.RandomRepeats, MinContigLength = o.MinContigLength };
        bool isdefDiff = o.MinCoverageDifference.HasValue;
        if (o.RandomRepeats)
        {
            cx.Rand = new GlibcRandom(o.RandomSeed ?? (uint)DateTimeOffset.UtcNow.ToUnixTimeSeconds());
            minCovDiff = 0;
            isdefDiff = true;
        }

        var file = MummerDeltaFile.Parse(deltaText, true);
        cx.Promer = file.DataType != "NUCMER";
        cx.MaxGap = o.MaxGap ?? (cx.Promer ? -1 : 1000);
        cx.MinCoverage = o.MinCoverage ?? (cx.Promer ? 50f : 95f);
        float defaultDiff = cx.Promer ? 30f : 10f;
        cx.MinCoverageDiff = isdefDiff ? minCovDiff : defaultDiff;
        cx.MinPidy = o.MinIdentity ?? (cx.Promer ? 55f : 90f);

        var contigs = ParseDelta(file, cx);
        TileContigs(contigs, cx);
        string? tigr = o.ContigList ? OutputContigs(contigs) : null;
        string? unusable = o.Unusable ? OutputUnusable(contigs) : null;
        LinkContigs(contigs, cx);
        string? pseudo = o.PseudoMolecule ? OutputPseudoMolecule(contigs, queryFasta ?? throw new ArgumentNullException(nameof(queryFasta))) : null;

        var sb = new StringBuilder();
        if (o.PrintAlignments)
            PrintTilingAlignments(sb, contigs);
        else if (o.Xml)
            PrintTilingXml(sb, contigs, file.QueryPath, o);
        else
            PrintTilingPath(sb, contigs);
        return new MummerShowTilingResult(sb.ToString(), string.Empty, tigr, unusable, pseudo);
    }

    private static List<QueryContig> ParseDelta(MummerDeltaFile file, Ctx cx)
    {
        var contigs = new List<QueryContig>();
        string currIdQ = string.Empty;
        foreach (var rec in file.Records)
        {
            if (!string.Equals(currIdQ, rec.QueryId, StringComparison.Ordinal))
            {
                currIdQ = rec.QueryId;
                contigs.Add(new QueryContig { IdQ = rec.QueryId, SeqLenQ = rec.QueryLength, TileLevel = rec.QueryLength < cx.MinContigLength ? IgnoreTileLevel : UnusedTileLevel });
            }

            var cur = contigs[^1];
            foreach (var a in rec.Alignments)
            {
                var st = new AlignStats { SeqLenR = rec.ReferenceLength, IdR = rec.ReferenceId, SR = a.StartR, ER = a.EndR, SQ = a.StartQ, EQ = a.EndQ };
                st.DirQ = (st.SR <= st.ER && st.SQ <= st.EQ) || (st.SR > st.ER && st.SQ > st.EQ) ? ForwardChar : ReverseChar;
                if (st.SR > st.ER)
                    (st.SR, st.ER) = (st.ER, st.SR);
                if (st.SQ > st.EQ)
                    (st.SQ, st.EQ) = (st.EQ, st.SQ);
                if (st.DirQ == ReverseChar)
                {
                    long temp = st.SQ;
                    st.SQ = RevC(st.EQ, rec.QueryLength);
                    st.EQ = RevC(temp, rec.QueryLength);
                }

                st.Idy = a.Identity;
                cur.Aligns.Add(st);
            }
        }

        foreach (var c in contigs)
            LibstdcxxSort.Sort(c.Aligns, (a, b) =>
            {
                int cmp = string.CompareOrdinal(a.IdR, b.IdR);
                if (cmp < 0)
                    return true;
                if (cmp > 0)
                    return false;
                return a.SQ < b.SQ;
            });
        return contigs;
    }

    private static void LongestConsistentSubset(List<AlignStats> list, int begin, int end, Ctx cx)
    {
        int n = end - begin;
        var score = new long[n];
        var from = new long[n];
        for (int i = 0; i < n; i++)
        {
            var ap = list[begin + i];
            if (ap.IsTiled)
                throw new InvalidOperationException("Assertion `Ap->isTiled == false' failed.");
            score[i] = ap.EQ - ap.SQ + 1;
            from[i] = -1;
        }

        for (int i = 0; i < n; i++)
        {
            var ai = list[begin + i];
            for (int j = 0; j < i; j++)
            {
                var aj = list[begin + j];
                if (ai.DirQ != aj.DirQ)
                    continue;
                long olapR = aj.ER - ai.SR + 1;
                if (cx.Circular && olapR > aj.SeqLenR / 2)
                    olapR = (aj.ER - aj.SeqLenR) + (1 - ai.SR);
                if (cx.MaxGap >= 0 && olapR < -cx.MaxGap)
                    continue;
                long olap = olapR < 0 ? 0 : olapR;
                long olapQ = aj.EQ - ai.SQ + 1;
                if (cx.MaxGap >= 0 && olapQ < -cx.MaxGap)
                    continue;
                olap = olap > olapQ ? olap : olapQ;
                if (score[j] + (ai.EQ - ai.SQ + 1) - olap > score[i])
                {
                    from[i] = j;
                    score[i] = score[j] + (ai.EQ - ai.SQ + 1) - olap;
                }
            }
        }

        int best = 0;
        for (int i = 1; i < n; i++)
        {
            if (score[i] > score[best])
                best = i;
            else if (score[i] == score[best] && list[begin + i].Idy > list[begin + best].Idy)
                best = i;
        }

        if (cx.RandomRepeats)
        {
            var ties = new List<int>();
            for (int i = 0; i < n; i++)
                if (score[i] == score[best] && (double)list[begin + i].Idy >= (double)list[begin + best].Idy - RepeatPidyDiff)
                    ties.Add(i);
            if (ties.Count > 1)
                best = ties[(int)((double)ties.Count * cx.Rand.Next() / (2147483647 + 1.0))];
        }

        for (long i = best; i >= 0; i = from[i])
            list[begin + (int)i].IsTiled = true;
    }

    private static void PlaceContig(QueryContig cp, Ctx cx)
    {
        AlignStats? ap = null;
        float maxCov = -float.MaxValue;
        foreach (var aip in cp.Aligns)
        {
            if (!aip.IsTiled)
                continue;
            if (ap == null)
                cp.LoTrim = aip.SQ - 1;
            cp.HiTrim = cp.SeqLenQ - aip.EQ;
            float cov = AlignmentQryCoverage(aip, cp.SeqLenQ);
            if (cov > maxCov)
            {
                maxCov = cov;
                ap = aip;
            }
        }

        if (ap == null)
            throw new InvalidOperationException("Assertion `Ap != Cp->Aligns.end( )' failed.");
        cp.DirQ = ap.DirQ;
        cp.IdR = ap.IdR;
        cp.SeqLenR = ap.SeqLenR;
        long start = ap.SQ;
        cp.StartR = ap.SR - (start - 1);
        cp.EndR = cp.StartR + cp.SeqLenQ - 1;
        if (cx.Circular && cp.EndR > cp.SeqLenR)
        {
            cp.StartR -= cp.SeqLenR;
            cp.EndR -= cp.SeqLenR;
        }

        if (cp.StartR > cp.EndR || cp.LoTrim < 0 || cp.HiTrim < 0)
            throw new InvalidOperationException("Assertion failed in placeContig.");
    }

    private static void TileContigs(List<QueryContig> contigs, Ctx cx)
    {
        foreach (var cp in contigs)
        {
            if (cp.TileLevel != UnusedTileLevel)
                continue;
            float maxCov = -float.MaxValue, maxxCov = -float.MaxValue, idy = -float.MaxValue;
            float maxCovHang = -float.MaxValue, maxxCovHang = -float.MaxValue, idyHang = -float.MaxValue;
            string idR = string.Empty, idRHang = string.Empty;
            var aligns = cp.Aligns;
            int ap = 0;
            while (ap < aligns.Count)
            {
                int aip;
                for (aip = ap + 1; aip < aligns.Count; aip++)
                    if (!string.Equals(aligns[ap].IdR, aligns[aip].IdR, StringComparison.Ordinal))
                        break;

                LongestConsistentSubset(aligns, ap, aip, cx);
                long hang = 0, end = -1, start = -1;
                for (int t = ap; t < aip; t++)
                {
                    if (aligns[t].IsTiled)
                    {
                        if (aligns[t].SQ - aligns[t].SR > 0)
                            hang += aligns[t].SQ - aligns[t].SR;
                        start = aligns[t].SR;
                        break;
                    }
                }

                for (int t = aip - 1; t >= ap; t--)
                {
                    if (aligns[t].IsTiled)
                    {
                        if (aligns[t].ER + (cp.SeqLenQ - aligns[t].EQ) > aligns[t].SeqLenR)
                            hang += (aligns[t].ER + (cp.SeqLenQ - aligns[t].EQ)) - aligns[t].SeqLenR;
                        end = aligns[t].ER;
                        break;
                    }
                }

                if (end <= 0 || start <= 0)
                    throw new InvalidOperationException("Assertion `end > 0 && start > 0' failed.");

                float cov = cx.Promer ? SubsetQrySyntenyCoverage(aligns, ap, aip, cp.SeqLenQ) : SubsetQryCoverage(aligns, ap, aip, cp.SeqLenQ);
                if (cov >= maxCov)
                {
                    float tmpidy = SubsetIdentity(aligns, ap, aip);
                    if (cov > maxCov || (cov == maxCov && tmpidy > idy))
                    {
                        idR = aligns[ap].IdR;
                        maxxCov = maxCov;
                        maxCov = cov;
                        idy = tmpidy;
                    }
                }
                else if (cov > maxxCov)
                {
                    maxxCov = cov;
                }

                if (!cx.Circular)
                {
                    float covHang = cx.Promer ? SubsetQrySyntenyCoverage(aligns, ap, aip, cp.SeqLenQ - hang) : SubsetQryCoverage(aligns, ap, aip, cp.SeqLenQ - hang);
                    if (covHang >= maxCovHang)
                    {
                        float tmpidy = SubsetIdentity(aligns, ap, aip);
                        if (covHang > maxCovHang || (covHang == maxCovHang && tmpidy > idyHang))
                        {
                            idRHang = aligns[ap].IdR;
                            maxxCovHang = maxCovHang;
                            maxCovHang = covHang;
                            idyHang = tmpidy;
                        }
                    }
                    else if (covHang > maxxCovHang)
                    {
                        maxxCovHang = covHang;
                    }
                }

                for (; ap < aip; ap++)
                {
                    var a = aligns[ap];
                    if (a.IsTiled)
                        continue;
                    if (start <= end && a.SR <= end && a.ER >= start)
                        continue;
                    if (start > end && (a.ER >= start || a.SR <= end))
                        continue;
                    cov = AlignmentQryCoverage(a, cp.SeqLenQ);
                    if (cov > maxxCov)
                        maxxCov = cov;
                    if (!cx.Circular)
                    {
                        hang = 0;
                        if (a.SQ - a.SR > 0)
                            hang += a.SQ - a.SR;
                        if (a.ER + (cp.SeqLenQ - a.EQ) > a.SeqLenR)
                            hang += (a.ER + (cp.SeqLenQ - a.EQ)) - a.SeqLenR;
                        float covHang = AlignmentQryCoverage(a, cp.SeqLenQ - hang);
                        if (covHang > maxxCovHang)
                            maxxCovHang = covHang;
                    }
                }
            }

            if (maxCov - maxxCov >= cx.MinCoverageDiff && maxCov >= cx.MinCoverage && idy >= cx.MinPidy)
            {
                foreach (var a in aligns)
                    if (a.IsTiled && !string.Equals(a.IdR, idR, StringComparison.Ordinal))
                        a.IsTiled = false;
                cp.TileLevel = UsedTileLevel;
                PlaceContig(cp, cx);
            }
            else if (!cx.Circular && maxCovHang - maxxCovHang >= cx.MinCoverageDiff && maxCovHang >= cx.MinCoverage && idyHang >= cx.MinPidy)
            {
                foreach (var a in aligns)
                    if (a.IsTiled && !string.Equals(a.IdR, idRHang, StringComparison.Ordinal))
                        a.IsTiled = false;
                cp.TileLevel = UsedTileLevel;
                PlaceContig(cp, cx);
            }
            else
            {
                cp.TileLevel = UnusableTileLevel;
            }
        }

        LibstdcxxSort.Sort(contigs, (a, b) =>
        {
            int cmp = string.CompareOrdinal(a.IdR, b.IdR);
            if (cmp < 0)
                return true;
            if (cmp > 0)
                return false;
            if (a.StartR < b.StartR)
                return true;
            if (a.StartR > b.StartR)
                return false;
            return a.SeqLenQ > b.SeqLenQ;
        });
    }

    private static void LinkContigs(List<QueryContig> contigs, Ctx cx)
    {
        int n = contigs.Count;
        int firstCp = n;
        int cp = -1;
        while (++cp < n)
        {
            if (contigs[cp].TileLevel != UsedTileLevel)
                continue;
            if (firstCp == n)
            {
                firstCp = cp;
                contigs[firstCp].IsLinkHead = true;
            }

            int nxtCp = n;
            long bestEnd = -long.MaxValue;
            for (int cip = cp + 1; cip < n; cip++)
            {
                if (contigs[cip].TileLevel != UsedTileLevel)
                    continue;
                if (!string.Equals(contigs[cip].IdR, contigs[cp].IdR, StringComparison.Ordinal))
                    break;
                if (contigs[cip].StartR > contigs[cp].EndR)
                {
                    if (nxtCp == n)
                        nxtCp = cip;
                    break;
                }

                contigs[cip].TileLevel = IgnoreTileLevel;
                if (contigs[cip].EndR > contigs[cp].EndR && contigs[cip].EndR > bestEnd)
                {
                    bestEnd = contigs[cip].EndR;
                    nxtCp = cip;
                }
            }

            if (nxtCp != n)
            {
                contigs[nxtCp].TileLevel = UsedTileLevel;
                contigs[cp].LinksTo = contigs[nxtCp];
            }
            else
            {
                if (cx.Circular)
                {
                    if (contigs[firstCp].StartR <= 0)
                    {
                        nxtCp = cp;
                        while (cp > firstCp)
                        {
                            if (contigs[cp].TileLevel == UsedTileLevel)
                            {
                                if (contigs[firstCp].SeqLenR + contigs[firstCp].StartR <= contigs[cp].EndR)
                                {
                                    nxtCp = cp;
                                    contigs[cp].TileLevel = IgnoreTileLevel;
                                }
                                else
                                {
                                    break;
                                }
                            }

                            cp--;
                        }

                        if (contigs[firstCp].SeqLenR + contigs[firstCp].StartR > contigs[nxtCp].StartR)
                        {
                            cp = nxtCp;
                            contigs[cp].TileLevel = UsedTileLevel;
                        }
                    }

                    contigs[cp].LinksTo = contigs[firstCp];
                }
                else
                {
                    contigs[cp].LinksTo = null;
                }

                firstCp = n;
            }
        }
    }

    private static string OutputContigs(List<QueryContig> contigsIn)
    {
        var contigs = new List<QueryContig>(contigsIn);
        LibstdcxxSort.Sort(contigs, (a, b) =>
        {
            int cmp = string.CompareOrdinal(a.IdR, b.IdR);
            if (cmp < 0)
                return true;
            if (cmp > 0)
                return false;
            return a.StartR + a.LoTrim < b.StartR + a.LoTrim;
        });
        var sb = new StringBuilder();
        var ci = CultureInfo.InvariantCulture;
        int cp = 0;
        while (cp < contigs.Count)
        {
            if (contigs[cp].TileLevel != UsedTileLevel)
            {
                cp++;
                continue;
            }

            long seqs = 0;
            int endCp;
            for (endCp = cp; endCp < contigs.Count; endCp++)
            {
                if (contigs[endCp].TileLevel != UsedTileLevel)
                    continue;
                if (!string.Equals(contigs[cp].IdR, contigs[endCp].IdR, StringComparison.Ordinal))
                    break;
                seqs++;
            }

            sb.Append("##").Append(contigs[cp].IdR).Append(' ').Append(seqs.ToString(ci)).Append(' ').Append(contigs[cp].SeqLenR.ToString(ci)).Append(" bases, 00000000 checksum.\n");
            for (; cp < endCp; cp++)
            {
                var c = contigs[cp];
                if (c.TileLevel != UsedTileLevel)
                    continue;
                if (c.DirQ == ForwardChar)
                    sb.Append('#').Append(c.IdQ).Append('(').Append((c.StartR + c.LoTrim - 1).ToString(ci)).Append(") [] ").Append(c.SeqLenQ.ToString(ci)).Append(" bases, 00000000 checksum. {")
                        .Append((c.LoTrim + 1).ToString(ci)).Append(' ').Append((c.SeqLenQ - c.HiTrim).ToString(ci)).Append("} <").Append((c.StartR + c.LoTrim).ToString(ci)).Append(' ').Append((c.EndR - c.HiTrim).ToString(ci)).Append(">\n");
                else
                    sb.Append('#').Append(c.IdQ).Append('(').Append((c.StartR + c.LoTrim - 1).ToString(ci)).Append(") [RC] ").Append(c.SeqLenQ.ToString(ci)).Append(" bases, 00000000 checksum. {")
                        .Append((c.SeqLenQ - c.LoTrim).ToString(ci)).Append(' ').Append((c.HiTrim + 1).ToString(ci)).Append("} <").Append((c.StartR + c.LoTrim).ToString(ci)).Append(' ').Append((c.EndR - c.HiTrim).ToString(ci)).Append(">\n");
            }
        }

        return sb.ToString();
    }

    private static string OutputUnusable(List<QueryContig> contigsIn)
    {
        var contigs = new List<QueryContig>(contigsIn);
        LibstdcxxSort.Sort(contigs, (a, b) => string.CompareOrdinal(a.IdQ, b.IdQ) < 0);
        var sb = new StringBuilder();
        foreach (var cp in contigs)
        {
            if (cp.TileLevel != UnusableTileLevel && cp.TileLevel != UnusedTileLevel)
                continue;
            foreach (var ap in cp.Aligns)
                PrintAlignment(sb, cp, ap);
        }

        return sb.ToString();
    }

    private static string OutputPseudoMolecule(List<QueryContig> contigs, string queryFasta)
    {
        foreach (var (id, seq) in NucmerPostnuc.ReadFasta(queryFasta, out _))
        {
            var cp = contigs.FirstOrDefault(c => c.TileLevel == UsedTileLevel && string.Equals(id, c.IdQ, StringComparison.Ordinal));
            if (cp == null)
                continue;
            if (seq.Length != cp.SeqLenQ)
                throw new InvalidOperationException("Assertion `strlen(A+1) == Cp->SeqLenQ' failed.");
            var bytes = new byte[seq.Length + 2];
            for (int i = 0; i < seq.Length; i++)
                bytes[i + 1] = (byte)seq[i];
            cp.SeqQ = cp.DirQ == ReverseChar ? MummerSequence.ReverseComplement(bytes) : bytes;
        }

        var sb = new StringBuilder();
        long ct = 0;
        void Put(char c)
        {
            sb.Append(c);
            if (++ct == 60)
            {
                ct = 0;
                sb.Append('\n');
            }
        }

        foreach (var begin in contigs)
        {
            if (!begin.IsLinkHead)
                continue;
            long start = 1;
            QueryContig? cp = begin;
            sb.Append(">pseudo_used_").Append(cp.IdR).Append('\n');
            while (cp != null)
            {
                long gap = 0;
                long end = cp.SeqLenQ;
                if (cp.LinksTo != null && !cp.LinksTo.IsLinkHead)
                {
                    long startR = cp.LinksTo.StartR;
                    long endR = cp.EndR;
                    gap = startR - endR - 1;
                    if (gap < 0)
                    {
                        int apb = 0;
                        var next = cp.LinksTo;
                        while (apb < next.Aligns.Count && !next.Aligns[apb].IsTiled)
                            apb++;
                        if (apb >= next.Aligns.Count)
                            throw new InvalidOperationException("Assertion `Apb < Cp->linksTo->Aligns.end( )' failed.");
                        int ape = cp.Aligns.Count - 1;
                        while (ape >= 0 && !cp.Aligns[ape].IsTiled)
                            ape--;
                        if (ape < 0)
                            throw new InvalidOperationException("Assertion `Ape < Cp->Aligns.rend( )' failed.");
                        if (cp.SeqLenQ - cp.Aligns[ape].EQ > next.Aligns[apb].SQ - 1
                            || (cp.SeqLenQ - cp.Aligns[ape].EQ == next.Aligns[apb].SQ - 1 && SubsetIdentity(cp.Aligns, 0, cp.Aligns.Count) < SubsetIdentity(next.Aligns, 0, next.Aligns.Count)))
                            end += gap;
                    }
                }

                if (cp.SeqQ == null)
                    throw new FormatException("ERROR: Sequence \"" + cp.IdQ + "\" was not found in the query file.");
                for (long i = start; i <= end; i++)
                    Put(char.ToUpperInvariant((char)cp.SeqQ[i]));
                for (long i = 1; i <= gap; i++)
                    Put('N');
                if (gap < 0 && end == cp.SeqLenQ)
                    start = 1 + -gap;
                cp = cp.LinksTo == null || cp.LinksTo.IsLinkHead ? null : cp.LinksTo;
            }

            if (ct != 0)
                sb.Append('\n');
            ct = 0;
        }

        return sb.ToString();
    }

    private static void PrintAlignment(StringBuilder sb, QueryContig cp, AlignStats ap)
    {
        var ci = CultureInfo.InvariantCulture;
        long len1 = ap.ER - ap.SR + 1, len2 = ap.EQ - ap.SQ + 1;
        float covA = Pct(len1, ap.SeqLenR);
        float covB = AlignmentQryCoverage(ap, cp.SeqLenQ);
        sb.Append(ap.SR.ToString(ci)).Append('\t').Append(ap.ER.ToString(ci)).Append('\t');
        if (ap.DirQ == ForwardChar)
            sb.Append(ap.SQ.ToString(ci)).Append('\t').Append(ap.EQ.ToString(ci)).Append('\t');
        else
            sb.Append(RevC(ap.SQ, cp.SeqLenQ).ToString(ci)).Append('\t').Append(RevC(ap.EQ, cp.SeqLenQ).ToString(ci)).Append('\t');
        sb.Append(len1.ToString(ci)).Append('\t').Append(len2.ToString(ci)).Append('\t').Append(CFormat.Fixed(ap.Idy, 2)).Append('\t');
        sb.Append(ap.SeqLenR.ToString(ci)).Append('\t').Append(cp.SeqLenQ.ToString(ci)).Append('\t').Append(CFormat.Fixed(covA, 2)).Append('\t').Append(CFormat.Fixed(covB, 2)).Append('\t');
        sb.Append(ap.IdR).Append('\t').Append(cp.IdQ).Append('\n');
    }

    private static void PrintTilingAlignments(StringBuilder sb, List<QueryContig> contigs)
    {
        foreach (var begin in contigs)
        {
            if (!begin.IsLinkHead)
                continue;
            QueryContig? cp = begin;
            while (cp != null)
            {
                foreach (var ap in cp.Aligns)
                    if (ap.IsTiled)
                        PrintAlignment(sb, cp, ap);
                cp = cp.LinksTo == null || cp.LinksTo.IsLinkHead ? null : cp.LinksTo;
            }
        }
    }

    private static void PrintTilingPath(StringBuilder sb, List<QueryContig> contigs)
    {
        var ci = CultureInfo.InvariantCulture;
        foreach (var begin in contigs)
        {
            if (!begin.IsLinkHead)
                continue;
            QueryContig? cp = begin;
            sb.Append('>').Append(cp.IdR).Append(' ').Append(cp.SeqLenR.ToString(ci)).Append(" bases\n");
            while (cp != null)
            {
                long len = cp.SeqLenQ;
                float pcov = SubsetQryCoverage(cp.Aligns, 0, cp.Aligns.Count, cp.SeqLenQ);
                float pidy = SubsetIdentity(cp.Aligns, 0, cp.Aligns.Count);
                long gap;
                if (cp.LinksTo == null)
                    gap = cp.SeqLenR - cp.EndR;
                else if (cp.LinksTo.IsLinkHead)
                    gap = (cp.SeqLenR + cp.LinksTo.StartR) - cp.EndR - 1;
                else
                    gap = cp.LinksTo.StartR - cp.EndR - 1;
                sb.Append(cp.StartR.ToString(ci)).Append('\t').Append(cp.EndR.ToString(ci)).Append('\t').Append(gap.ToString(ci)).Append('\t').Append(len.ToString(ci)).Append('\t')
                    .Append(CFormat.Fixed(pcov, 2)).Append('\t').Append(CFormat.Fixed(pidy, 2)).Append('\t').Append(cp.DirQ).Append('\t').Append(cp.IdQ).Append('\n');
                cp = cp.LinksTo == null || cp.LinksTo.IsLinkHead ? null : cp.LinksTo;
            }
        }
    }

    private static void PrintTilingXml(StringBuilder sb, List<QueryContig> contigs, string queryName, MummerShowTilingOptions o)
    {
        var ci = CultureInfo.InvariantCulture;
        DateTime now = o.Now ?? DateTime.Now;
        string timeStr = now.ToString("ddd MMM ", ci) + now.Day.ToString(ci).PadLeft(2) + now.ToString(" HH:mm:ss yyyy", ci);
        sb.Append("<?xml version = \"1.0\" ?>\n<EVIDENCE ID = \"project_1\"\n        DATE = \"").Append(timeStr).Append("\"\n     PROJECT = \"").Append(queryName).Append("\"\n  PARAMETERS = \"")
            .Append(string.Join(' ', o.CommandLine ?? [])).Append("\">\n");
        foreach (var begin in contigs)
        {
            if (!begin.IsLinkHead)
                continue;
            QueryContig? cp = begin;
            while (cp != null)
            {
                sb.Append("\t<CONTIG ID = \"contig_").Append(cp.IdQ).Append("\" NAME = \"").Append(cp.IdQ).Append("\" LEN = \"").Append(cp.SeqLenQ.ToString(ci)).Append("\"/>\n");
                cp = cp.LinksTo == null || cp.LinksTo.IsLinkHead ? null : cp.LinksTo;
            }
        }

        sb.Append('\n');
        long ct = 0;
        foreach (var begin in contigs)
        {
            if (!begin.IsLinkHead)
                continue;
            QueryContig? cp = begin;
            while (cp != null)
            {
                if (cp.LinksTo == null)
                    break;
                long gap = cp.LinksTo.IsLinkHead ? (cp.SeqLenR + cp.LinksTo.StartR) - cp.EndR - 1 : cp.LinksTo.StartR - cp.EndR - 1;
                sb.Append("\t<LINK ID = \"link_").Append((++ct).ToString(ci)).Append("\" SIZE = \"").Append(gap.ToString(ci)).Append("\" TYPE = \"MUMmer\">\n");
                sb.Append("\t\t<CONTIG ID = \"contig_").Append(cp.IdQ).Append("\" ORI = \"").Append(cp.DirQ == ForwardChar ? "BE" : "EB").Append("\">\n");
                foreach (var ap in cp.Aligns)
                {
                    if (!ap.IsTiled)
                        continue;
                    sb.Append("\t\t");
                    PrintAlignment(sb, cp, ap);
                }

                sb.Append("\t\t</CONTIG>\n");
                sb.Append("\t\t<CONTIG ID = \"contig_").Append(cp.LinksTo.IdQ).Append("\" ORI = \"").Append(cp.LinksTo.DirQ == ForwardChar ? "BE" : "EB").Append("\">\n");
                foreach (var ap in cp.LinksTo.Aligns)
                {
                    if (!ap.IsTiled)
                        continue;
                    sb.Append("\t\t");
                    PrintAlignment(sb, cp.LinksTo, ap);
                }

                sb.Append("\t\t</CONTIG>\n\t</LINK>\n");
                cp = cp.LinksTo.IsLinkHead ? null : cp.LinksTo;
            }
        }

        sb.Append("</EVIDENCE>\n");
    }
}
