using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>Options of <see cref="MummerShowCoords.Run"/> (the switches of <c>show-coords</c> 3.23).</summary>
public sealed record MummerShowCoordsOptions
{
    /// <summary><c>-b</c>: merge overlapping alignments regardless of direction or frame, no identity columns.</summary>
    public bool Brief { get; init; }

    /// <summary><c>-B</c>: btab output.</summary>
    public bool Btab { get; init; }

    /// <summary><c>-c</c>: percent coverage columns.</summary>
    public bool ShowCoverage { get; init; }

    /// <summary><c>-d</c>: the alignment direction in the FRM columns (always for promer data).</summary>
    public bool ShowDirection { get; init; }

    /// <summary><c>-H</c> inverted: print the header (default true).</summary>
    public bool PrintHeader { get; init; } = true;

    /// <summary><c>-I</c>: minimum alignment identity (0 - 100).</summary>
    public float MinIdentity { get; init; }

    /// <summary><c>-k</c>: knockout (promer data only).</summary>
    public bool Knockout { get; init; }

    /// <summary><c>-l</c>: sequence length columns.</summary>
    public bool ShowSeqLengths { get; init; }

    /// <summary><c>-L</c>: minimum alignment length.</summary>
    public long MinLength { get; init; }

    /// <summary><c>-o</c>: annotate overlaps ([IDENTITY] [CONTAINED] [CONTAINS] [BEGIN] [END]).</summary>
    public bool AnnotateOverlaps { get; init; }

    /// <summary><c>-q</c>: sort by query.</summary>
    public bool SortByQuery { get; init; }

    /// <summary><c>-r</c>: sort by reference.</summary>
    public bool SortByReference { get; init; }

    /// <summary><c>-g</c> (deprecated): longest ascending subset.</summary>
    public bool Las { get; init; }

    /// <summary><c>-G</c> (hidden): weighted longest ascending subset.</summary>
    public bool WeightedLas { get; init; }

    /// <summary><c>-T</c>: tab-delimited output.</summary>
    public bool Tabular { get; init; }

    /// <summary><c>-w</c> (hidden): warnings column ([DUPLICATE] [CONTAINS] [SHADOWED] [OVERLAPS]).</summary>
    public bool ShowWarnings { get; init; }

    /// <summary>The <c>srand (time (NULL))</c> seed of the program (it breaks ties of <c>-g</c> / <c>-G</c> with <c>rand ()</c>); null = the current time.</summary>
    public uint? RandomSeed { get; init; }

    /// <summary>The date a btab line shows (<c>%b %d %Y</c> of the current time); null = now.</summary>
    public DateTime? Now { get; init; }
}

/// <summary>Output of a MUMmer text program: standard output and standard error.</summary>
/// <param name="Text">What the program writes to standard output.</param>
/// <param name="Messages">What it writes to standard error (warnings).</param>
public sealed record MummerTextOutput(string Text, string Messages);

/// <summary>
/// <c>show-coords</c> of MUMmer 3.23 (<c>src/tigr/show-coords.cc</c>): the coordinates of the alignments of a delta file in the human-readable table, the tab-delimited table (-T) or btab (-B), with
/// the identity / coverage / length / frame columns, filters (-I -L), sorting (-r -q), merging (-b), overlap annotation (-o) and the longest-ascending-subset options.
/// </summary>
public static class MummerShowCoords
{
    private sealed class AlignStats
    {
        public bool IsLas;
        public float Idy, Sim, Stp;
        public int FrameA, FrameB;
        public long SA, EA, SB, EB, SeqLenA, SeqLenB;
        public string IdA = string.Empty, IdB = string.Empty;
        public string Annot = string.Empty;
    }

    private sealed class LasStats
    {
        public bool Ori;
        public long Score, WScore, From, WFrom;
        public bool WPoint;
        public AlignStats Asi = null!;
    }

    private const float MinOverlapP = 0.5f;
    private const float MinLengthP = 0.75f;
    private const float MaxOhangP = 0.05f;
    private const float MaxPidyDiff = 0.01f;

    private static long RevC(long coord, long len)
    {
        if (len - coord + 1 <= 0)
            throw new InvalidOperationException("Assertion `Len - Coord + 1 > 0' failed.");
        return len - coord + 1;
    }

    /// <summary>Runs <c>show-coords</c> on the text of a delta file.</summary>
    /// <param name="deltaText">The delta file.</param>
    /// <param name="options">The switches (null = defaults).</param>
    /// <exception cref="ArgumentOutOfRangeException">A percentage outside [0, 100] or a negative length (the program's exits).</exception>
    /// <exception cref="FormatException">The delta file cannot be parsed.</exception>
    public static MummerTextOutput Run(string deltaText, MummerShowCoordsOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(deltaText);
        var o = options ?? new MummerShowCoordsOptions();
        var err = new StringBuilder();
        bool isBrief = o.Brief, isBtab = o.Btab, isShowCoverage = o.ShowCoverage, isShowDir = o.ShowDirection, isPrintHeader = o.PrintHeader;
        bool isKnockout = o.Knockout, isShowSeqLens = o.ShowSeqLengths, isAnnotateOverlaps = o.AnnotateOverlaps, isSortByQuery = o.SortByQuery, isSortByReference = o.SortByReference;
        bool isLas = o.Las || o.WeightedLas, isWLas = o.WeightedLas, isTabular = o.Tabular, isShowWarnings = o.ShowWarnings;

        if (isShowWarnings && isAnnotateOverlaps)
        {
            err.Append("WARNING: -o and -w are mutually exclusive, -w ignored\n");
            isShowWarnings = false;
        }

        if (isWLas || isShowWarnings)
            err.Append("WARNING: hidden option used, I hope you know what you are doing\n");
        if (isLas && !isSortByReference && !isSortByQuery)
        {
            err.Append("WARNING: -g option requires -r or -q, using -r\n");
            isSortByReference = true;
        }

        if (isBtab && (!isPrintHeader || isShowCoverage || isShowSeqLens || isTabular || isShowWarnings || isBrief || isAnnotateOverlaps))
        {
            err.Append("WARNING: the options -bdHcloTw are all ignored for btab\n");
            isPrintHeader = isShowCoverage = isShowSeqLens = isTabular = isShowWarnings = isBrief = isAnnotateOverlaps = false;
        }

        if (isBrief && (isKnockout || isLas))
        {
            err.Append("WARNING: -k and -g have no effect when used with -b, -k and -g ignored\n");
            isKnockout = false;
            isLas = false;
        }

        if (isSortByReference && isSortByQuery)
            err.Append("WARNING: both -r and -q were passed, -q ignored\n");
        if (o.MinIdentity < 0.0 || o.MinIdentity > 100.0)
            throw new ArgumentOutOfRangeException(nameof(options), "ERROR: Percentages must be within the range [0.0. 100.0]");
        if (o.MinLength < 0)
            throw new ArgumentOutOfRangeException(nameof(options), "ERROR: Length values must be positive");
        bool isAnnotation = isShowWarnings || isAnnotateOverlaps;

        var file = MummerDeltaFile.Parse(deltaText, true);
        bool promer = file.DataType != "NUCMER";
        var stats = new List<AlignStats>();
        foreach (var rec in file.Records)
        {
            foreach (var a in rec.Alignments)
            {
                var cur = new AlignStats
                {
                    SeqLenA = rec.ReferenceLength, SeqLenB = rec.QueryLength, IdA = rec.ReferenceId, IdB = rec.QueryId,
                    SA = a.StartR, EA = a.EndR, SB = a.StartQ, EB = a.EndQ, Idy = a.Identity, Sim = a.Similarity, Stp = a.Stop,
                };
                if (cur.Idy < o.MinIdentity)
                    continue;
                if ((double)cur.Idy > 99.99 && cur.Idy != 100)
                    cur.Idy = (float)99.99;
                if ((double)cur.Sim > 99.99 && cur.Sim != 100)
                    cur.Sim = (float)99.99;
                if ((double)cur.Stp > 99.99 && cur.Stp != 100)
                    cur.Stp = (float)99.99;
                if ((double)cur.Stp < 0.01 && cur.Stp != 0)
                    cur.Stp = (float)0.01;
                if (Math.Abs(cur.SA - cur.EA) + 1 < o.MinLength || Math.Abs(cur.SB - cur.EB) + 1 < o.MinLength)
                    continue;

                long sA = cur.SA, eA = cur.EA, sB = cur.SB, eB = cur.EB;
                int frameA = 1, frameB = 1;
                if (sA > eA)
                {
                    sA = RevC(sA, cur.SeqLenA);
                    _ = RevC(eA, cur.SeqLenA);
                    frameA += 3;
                }

                if (sB > eB)
                {
                    sB = RevC(sB, cur.SeqLenB);
                    _ = RevC(eB, cur.SeqLenB);
                    frameB += 3;
                }

                if (isBrief)
                {
                    cur.FrameA = frameA > 3 ? (frameA - 3) * -1 : frameA;
                    cur.FrameB = frameB > 3 ? (frameB - 3) * -1 : frameB;
                    stats.Add(cur);
                    continue;
                }

                if (!promer)
                {
                    if (frameA != 1 || (frameB != 1 && frameB != 4))
                        throw new InvalidOperationException("Assertion `frameA == 1 && (frameB == 1 || frameB == 4)' failed.");
                }
                else
                {
                    frameA += (int)((sA + 2) % 3);
                    frameB += (int)((sB + 2) % 3);
                }

                cur.FrameA = frameA > 3 ? (frameA - 3) * -1 : frameA;
                cur.FrameB = frameB > 3 ? (frameB - 3) * -1 : frameB;
                stats.Add(cur);
            }
        }

        if (isKnockout && !promer)
        {
            err.Append("WARNING: -k only applies to PROmer data, -k ignored\n");
            isKnockout = false;
        }

        if (promer)
            isShowDir = true;
        if (isBrief || isKnockout)
            SimplifyAlignments(stats, isBrief, isKnockout);
        if (isLas)
            FlagLas(stats, isWLas, new GlibcRandom(o.RandomSeed ?? (uint)DateTimeOffset.UtcNow.ToUnixTimeSeconds()));

        if (isSortByReference)
            LibstdcxxSort.Sort(stats, (a, b) => IdASAIdBSBLess(a, b));
        else if (isSortByQuery)
            LibstdcxxSort.Sort(stats, (a, b) => IdBSBIdASALess(a, b));

        if (isShowWarnings)
            GenerateWarnings(stats);
        if (isAnnotateOverlaps)
            AnnotateOverlapsOf(stats);

        var sb = new StringBuilder();
        if (isBtab)
            PrintBtab(sb, stats, file, promer, o.Now ?? DateTime.Now);
        else if (isTabular)
            PrintTabular(sb, stats, file, promer, isBrief, isShowSeqLens, isShowCoverage, isShowDir, isPrintHeader, isAnnotation);
        else
            PrintHuman(sb, stats, file, promer, isBrief, isShowSeqLens, isShowCoverage, isShowDir, isPrintHeader, isAnnotation);
        return new MummerTextOutput(sb.ToString(), err.ToString());
    }

    private static bool IdASAIdBSBLess(AlignStats a, AlignStats b)
    {
        int i = string.CompareOrdinal(a.IdA, b.IdA);
        int j = string.CompareOrdinal(a.IdB, b.IdB);
        if (i < 0)
            return true;
        if (i > 0)
            return false;
        if (a.SA < b.SA)
            return true;
        if (a.SA > b.SA)
            return false;
        if (j < 0)
            return true;
        if (j > 0)
            return false;
        return a.SB < b.SB;
    }

    private static bool IdBSBIdASALess(AlignStats a, AlignStats b)
    {
        int i = string.CompareOrdinal(a.IdB, b.IdB);
        int j = string.CompareOrdinal(a.IdA, b.IdA);
        if (i < 0)
            return true;
        if (i > 0)
            return false;
        if (a.SB < b.SB)
            return true;
        if (a.SB > b.SB)
            return false;
        if (j < 0)
            return true;
        if (j > 0)
            return false;
        return a.SA < b.SA;
    }

    private static void Swap(ref long a, ref long b) => (a, b) = (b, a);

    private static void SimplifyAlignments(List<AlignStats> stats, bool isBrief, bool isKnockout)
    {
        int sp = 0;
        while (sp < stats.Count)
        {
            int sep;
            for (sep = sp; sep < stats.Count; sep++)
            {
                if (stats[sp].IdA != stats[sep].IdA || stats[sp].IdB != stats[sep].IdB)
                    break;
                var s = stats[sep];
                if (s.FrameA < 0)
                    Swap(ref s.SA, ref s.EA);
                if (s.FrameB < 0)
                    Swap(ref s.SB, ref s.EB);
                if (isBrief)
                {
                    s.FrameA = 1;
                    s.FrameB = 1;
                }
            }

            var range = stats.GetRange(sp, sep - sp);
            LibstdcxxSort.Sort(range, (a, b) => a.SA < b.SA);
            for (int k = 0; k < range.Count; k++)
                stats[sp + k] = range[k];

            for (; sp < sep; sp++)
            {
                var p = stats[sp];
                if (p.SA < 0)
                    continue;
                for (int sip = sp + 1; sip < sep; sip++)
                {
                    var q = stats[sip];
                    if (q.SA > p.EA)
                        break;
                    if (isBrief)
                    {
                        if (q.SB <= p.EB && q.EB >= p.SB)
                        {
                            p.SA = Math.Min(p.SA, q.SA);
                            p.SB = Math.Min(p.SB, q.SB);
                            p.EA = Math.Max(p.EA, q.EA);
                            p.EB = Math.Max(p.EB, q.EB);
                            q.SA *= -1;
                        }
                    }
                    else if (isKnockout)
                    {
                        if ((long)(p.FrameA * p.FrameB) * (q.FrameA * q.FrameB) < 0)
                            continue;
                        long lenA = p.EA - p.SA + 1, lenAi = q.EA - q.SA + 1, lenB = p.EB - p.SB + 1, lenBi = q.EB - q.SB + 1;
                        long olapA = p.EA - q.SA + 1;
                        long olapB = p.EB - p.SB + 1;
                        long tmp = p.EB - q.EB;
                        if (tmp > 0)
                            olapB -= tmp;
                        tmp = q.SB - p.SB;
                        if (tmp > 0)
                            olapB -= tmp;
                        float fa = (float)olapA / (float)(lenA < lenAi ? lenA : lenAi);
                        float fb = (float)olapB / (float)(lenB < lenBi ? lenB : lenBi);
                        if (fa >= MinOverlapP && fb >= MinOverlapP)
                        {
                            float la = (float)lenA / (float)lenAi, lb = (float)lenB / (float)lenBi;
                            if (p.Sim > q.Sim && la >= MinLengthP && lb >= MinLengthP)
                            {
                                q.SA *= -1;
                            }
                            else
                            {
                                p.SA *= -1;
                                break;
                            }
                        }
                    }
                    else
                    {
                        throw new InvalidOperationException("Assertion `false' failed.");
                    }
                }

                if (p.SA > 0 && p.FrameA < 0)
                    Swap(ref p.SA, ref p.EA);
                if (p.SA > 0 && p.FrameB < 0)
                    Swap(ref p.SB, ref p.EB);
            }
        }
    }

    private static void FlagLas(List<AlignStats> stats, bool isWLas, GlibcRandom rand)
    {
        int sp = 0;
        while (sp < stats.Count)
        {
            var las = new List<LasStats>();
            int sep;
            for (sep = sp; sep < stats.Count; sep++)
            {
                var s = stats[sep];
                if (stats[sp].IdA != s.IdA || stats[sp].IdB != s.IdB)
                    break;
                s.IsLas = false;
                if (s.SA < 0)
                    continue;
                if (s.FrameA < 0)
                    Swap(ref s.SA, ref s.EA);
                if (s.FrameB < 0)
                    Swap(ref s.SB, ref s.EB);
                var alas = new LasStats();
                if ((s.FrameA < 0 && s.FrameB < 0) || (s.FrameA > 0 && s.FrameB > 0))
                {
                    alas.Ori = false;
                }
                else
                {
                    alas.Ori = true;
                    Swap(ref s.SB, ref s.EB);
                    s.SB = RevC(s.SB, s.SeqLenB);
                    s.EB = RevC(s.EB, s.SeqLenB);
                }

                long lenA = s.EA - s.SA + 1, lenB = s.EB - s.SB + 1;
                alas.Score = alas.WScore = lenA < lenB ? lenA : lenB;
                alas.From = alas.WFrom = -1;
                alas.WPoint = false;
                alas.Asi = s;
                las.Add(alas);
            }

            LibstdcxxSort.Sort(las, (a, b) => a.Asi.SA < b.Asi.SA);
            int n = las.Count;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    if (las[i].Ori != las[j].Ori)
                        continue;
                    long lenA = las[i].Asi.EA - las[i].Asi.SA + 1, lenB = las[i].Asi.EB - las[i].Asi.SB + 1;
                    long len = lenA < lenB ? lenA : lenB;
                    long olapA = las[j].Asi.EA - las[i].Asi.SA + 1;
                    long olap = olapA > 0 ? olapA : 0;
                    if (las[j].Score + len - olap > las[i].WScore)
                    {
                        las[i].WFrom = j;
                        las[i].WScore = las[j].Score + len - olap;
                        las[i].WPoint = true;
                    }

                    long olapB = las[j].Asi.EB - las[i].Asi.SB + 1;
                    olap = olapB > olap ? olapB : olap;
                    if (las[j].Score + len - olap > las[i].Score)
                    {
                        las[i].From = j;
                        las[i].Score = las[j].Score + len - olap;
                    }

                    if (las[j].WScore + len - olap > las[i].WScore)
                    {
                        las[i].WFrom = j;
                        las[i].WScore = las[j].WScore + len - olap;
                        las[i].WPoint = false;
                    }
                }
            }

            int best = 0;
            for (int i = 0; i < n; i++)
            {
                if (isWLas
                    ? las[i].WScore > las[best].WScore || (las[i].WScore == las[best].WScore && las[i].Asi.Idy > las[best].Asi.Idy)
                    : las[i].Score > las[best].Score || (las[i].Score == las[best].Score && las[i].Asi.Idy > las[best].Asi.Idy))
                    best = i;
            }

            var ties = new List<int>();
            for (int i = 0; i < n; i++)
            {
                var a = las[i].Asi;
                if (las[i].Ori)
                {
                    Swap(ref a.SB, ref a.EB);
                    a.SB = RevC(a.SB, a.SeqLenB);
                    a.EB = RevC(a.EB, a.SeqLenB);
                }

                if (a.FrameA < 0)
                    Swap(ref a.SA, ref a.EA);
                if (a.FrameB < 0)
                    Swap(ref a.SB, ref a.EB);
                if (((isWLas && las[i].WScore == las[best].WScore) || (!isWLas && las[i].Score == las[best].Score))
                    && (double)a.Idy >= (double)las[best].Asi.Idy - MaxPidyDiff)
                    ties.Add(i);
            }

            if (ties.Count > 1)
                best = ties[(int)((double)ties.Count * rand.Next() / (2147483647 + 1.0))];

            if (isWLas)
            {
                bool wrap = false;
                long next;
                for (long i = best; i >= 0; i = next)
                {
                    las[(int)i].Asi.IsLas = true;
                    if (wrap)
                    {
                        next = las[(int)i].From;
                        las[(int)i].WPoint = false;
                    }
                    else
                    {
                        next = las[(int)i].WFrom;
                    }

                    if (las[(int)i].WPoint)
                        wrap = true;
                }
            }
            else
            {
                for (long i = best; i >= 0; i = las[(int)i].From)
                    las[(int)i].Asi.IsLas = true;
            }

            sp = sep;
        }

        foreach (var s in stats)
            if (s.SA > 0 && !s.IsLas)
                s.SA *= -1;
    }

    private static void AnnotateOverlapsOf(List<AlignStats> stats)
    {
        foreach (var s in stats)
        {
            s.Annot = string.Empty;
            if (s.SA < 0)
                continue;
            long loR = s.SA, hiR = s.EA, loQ = s.SB, hiQ = s.EB;
            if (s.FrameA < 0)
                Swap(ref loR, ref hiR);
            if (s.FrameB < 0)
                Swap(ref loQ, ref hiQ);
            if (!(s.FrameA < 0 && s.FrameB < 0) && !(s.FrameA > 0 && s.FrameB > 0))
            {
                Swap(ref loQ, ref hiQ);
                loQ = RevC(loQ, s.SeqLenB);
                hiQ = RevC(hiQ, s.SeqLenB);
            }

            long span = hiR - loR > hiQ - loQ ? hiR - loR + 1 : hiQ - loQ + 1;
            float scaled = MaxOhangP * (float)span;
            long maxOhang = (long)scaled;
            if ((loR <= maxOhang || loQ <= maxOhang) && (hiR + maxOhang > s.SeqLenA || hiQ + maxOhang > s.SeqLenB))
            {
                if (loR <= maxOhang && loQ <= maxOhang && hiR + maxOhang > s.SeqLenA && hiQ + maxOhang > s.SeqLenB)
                    s.Annot = "[IDENTITY]";
                else if (loR <= maxOhang && hiR + maxOhang > s.SeqLenA)
                    s.Annot = "[CONTAINED]";
                else if (loQ <= maxOhang && hiQ + maxOhang > s.SeqLenB)
                    s.Annot = "[CONTAINS]";
                else if (loQ <= maxOhang)
                    s.Annot = "[END]";
                else if (loR <= maxOhang)
                    s.Annot = "[BEGIN]";
                else
                    throw new InvalidOperationException("Assertion `false' failed.");
            }
        }
    }

    private static void GenerateWarnings(List<AlignStats> stats)
    {
        if (stats.Count > 0)
            stats[0].Annot = string.Empty;
        for (int sip = 1; sip < stats.Count; sip++)
        {
            var cur = stats[sip];
            cur.Annot = string.Empty;
            if (cur.SA < 0)
                continue;
            int sprev = sip - 1;
            while (sprev >= 0 && stats[sprev].SA < 0)
                sprev--;

            if (sprev < 0 || cur.IdA != stats[sprev].IdA || cur.IdB != stats[sprev].IdB)
                continue;
            var prev = stats[sprev];
            long ps1 = prev.SA, pe1 = prev.EA, ps2 = prev.SB, pe2 = prev.EB, s1 = cur.SA, e1 = cur.EA, s2 = cur.SB, e2 = cur.EB;
            if (prev.FrameA < 0)
                Swap(ref ps1, ref pe1);
            if (prev.FrameB < 0)
                Swap(ref ps2, ref pe2);
            if (cur.FrameA < 0)
                Swap(ref s1, ref e1);
            if (cur.FrameB < 0)
                Swap(ref s2, ref e2);

            if (s1 == ps1 && s2 == ps2 && e1 == pe1 && e2 == pe2)
                cur.Annot = "[DUPLICATE]";
            else if (s1 < ps1 && s2 < ps2 && e1 > pe1 && e2 > pe2)
                cur.Annot = "[CONTAINS]";
            else if (ps1 < s1 && ps2 < s2 && pe1 > e1 && pe2 > e2)
                cur.Annot = "[SHADOWED]";
            else if ((s1 >= ps1 && s1 <= pe1) || (e1 >= ps1 && e1 <= pe1))
                cur.Annot = "[OVERLAPS]";
            else if ((s2 >= ps2 && s2 <= pe2) || (e2 >= ps2 && e2 <= pe2))
                cur.Annot = "[OVERLAPS]";
        }
    }

    private static string L(long v, int width) => v.ToString(CultureInfo.InvariantCulture).PadLeft(width);

    private static string F2(double v, int width) => CFormat.Fixed(v, 2).PadLeft(width);

    private static float Coverage(long len, long seqLen)
    {
        float q = (float)len / (float)seqLen;
        return (float)((double)q * 100.0);
    }

    private static void PrintBtab(StringBuilder sb, List<AlignStats> stats, MummerDeltaFile file, bool promer, DateTime now)
    {
        string date = now.ToString("MMM dd yyyy", CultureInfo.InvariantCulture);
        string type = promer ? "PROMER" : "NUCMER";
        foreach (var s in stats)
        {
            if (s.SA < 0)
                continue;
            long len = Math.Abs(s.EB - s.SB) + 1;
            sb.Append(s.IdB).Append('\t').Append(date).Append('\t').Append(L(s.SeqLenB, 0)).Append('\t').Append(type).Append('\t').Append(file.ReferencePath).Append('\t').Append(s.IdA).Append('\t');
            sb.Append(L(s.SB, 0)).Append('\t').Append(L(s.EB, 0)).Append('\t').Append(L(s.SA, 0)).Append('\t').Append(L(s.EA, 0)).Append('\t')
                .Append(CFormat.Fixed(s.Idy, 6)).Append('\t').Append(CFormat.Fixed(s.Sim, 6)).Append('\t').Append(L(len, 0)).Append("\t0\t0\tNULL\t");
            sb.Append(promer ? s.FrameA : 0).Append('\t').Append(s.FrameB < 0 ? "Minus" : "Plus").Append('\t').Append(L(s.SeqLenA, 0)).Append("\t0\t0\n");
        }
    }

    private static void PrintHuman(StringBuilder sb, List<AlignStats> stats, MummerDeltaFile file, bool promer, bool isBrief, bool showLens, bool showCov, bool showDir, bool header, bool annotation)
    {
        if (header)
        {
            sb.Append(file.ReferencePath).Append(' ').Append(file.QueryPath).Append('\n').Append(promer ? "PROMER" : "NUCMER").Append("\n\n");
            sb.Append("[S1]".PadLeft(8)).Append(' ').Append("[E1]".PadLeft(8)).Append("  | ");
            sb.Append("[S2]".PadLeft(8)).Append(' ').Append("[E2]".PadLeft(8)).Append("  | ");
            sb.Append("[LEN 1]".PadLeft(8)).Append(' ').Append("[LEN 2]".PadLeft(8)).Append("  | ");
            if (!isBrief)
            {
                sb.Append("[% IDY]".PadLeft(8)).Append(' ');
                if (promer)
                    sb.Append("[% SIM]".PadLeft(8)).Append(' ').Append("[% STP]".PadLeft(8)).Append(' ');
                sb.Append(" | ");
            }

            if (showLens)
                sb.Append("[LEN R]".PadLeft(8)).Append(' ').Append("[LEN Q]".PadLeft(8)).Append("  | ");
            if (showCov)
                sb.Append("[COV R]".PadLeft(8)).Append(' ').Append("[COV Q]".PadLeft(8)).Append("  | ");
            if (!isBrief && showDir)
                sb.Append("[FRM]".PadLeft(5)).Append("  ");
            sb.Append("[TAGS]").Append('\n');
            if (showLens)
                sb.Append("=====================");
            if (showCov)
                sb.Append("=====================");
            if (!isBrief)
            {
                sb.Append("============");
                if (showDir)
                    sb.Append("=======");
                if (promer)
                    sb.Append("==================");
            }

            sb.Append("====================================").Append("=====================================\n");
        }

        foreach (var s in stats)
        {
            if (s.SA < 0)
                continue;
            long len1 = Math.Abs(s.EA - s.SA) + 1, len2 = Math.Abs(s.EB - s.SB) + 1;
            sb.Append(L(s.SA, 8)).Append(' ').Append(L(s.EA, 8)).Append("  | ").Append(L(s.SB, 8)).Append(' ').Append(L(s.EB, 8)).Append("  | ").Append(L(len1, 8)).Append(' ').Append(L(len2, 8)).Append("  | ");
            if (!isBrief)
            {
                sb.Append(F2(s.Idy, 8)).Append(' ');
                if (promer)
                    sb.Append(F2(s.Sim, 8)).Append(' ').Append(F2(s.Stp, 8)).Append(' ');
                sb.Append(" | ");
            }

            if (showLens)
                sb.Append(L(s.SeqLenA, 8)).Append(' ').Append(L(s.SeqLenB, 8)).Append("  | ");
            if (showCov)
                sb.Append(F2(Coverage(len1, s.SeqLenA), 8)).Append(' ').Append(F2(Coverage(len2, s.SeqLenB), 8)).Append("  | ");
            if (!isBrief && showDir)
                sb.Append(s.FrameA.ToString(CultureInfo.InvariantCulture).PadLeft(2)).Append(' ').Append(s.FrameB.ToString(CultureInfo.InvariantCulture).PadLeft(2)).Append("  ");
            sb.Append(s.IdA).Append('\t').Append(s.IdB);
            if (annotation)
                sb.Append('\t').Append(s.Annot);
            sb.Append('\n');
        }
    }

    private static void PrintTabular(StringBuilder sb, List<AlignStats> stats, MummerDeltaFile file, bool promer, bool isBrief, bool showLens, bool showCov, bool showDir, bool header, bool annotation)
    {
        if (header)
        {
            sb.Append(file.ReferencePath).Append(' ').Append(file.QueryPath).Append('\n').Append(promer ? "PROMER" : "NUCMER").Append("\n\n");
            sb.Append("[S1]\t[E1]\t[S2]\t[E2]\t[LEN 1]\t[LEN 2]\t");
            if (!isBrief)
            {
                sb.Append("[% IDY]\t");
                if (promer)
                    sb.Append("[% SIM]\t[% STP]\t");
            }

            if (showLens)
                sb.Append("[LEN R]\t[LEN Q]\t");
            if (showCov)
                sb.Append("[COV R]\t[COV Q]\t");
            if (!isBrief && showDir)
                sb.Append("[FRM]\t");
            sb.Append("[TAGS]\n");
        }

        foreach (var s in stats)
        {
            if (s.SA < 0)
                continue;
            long len1 = Math.Abs(s.EA - s.SA) + 1, len2 = Math.Abs(s.EB - s.SB) + 1;
            sb.Append(L(s.SA, 0)).Append('\t').Append(L(s.EA, 0)).Append('\t').Append(L(s.SB, 0)).Append('\t').Append(L(s.EB, 0)).Append('\t').Append(L(len1, 0)).Append('\t').Append(L(len2, 0)).Append('\t');
            if (!isBrief)
            {
                sb.Append(CFormat.Fixed(s.Idy, 2)).Append('\t');
                if (promer)
                    sb.Append(CFormat.Fixed(s.Sim, 2)).Append('\t').Append(CFormat.Fixed(s.Stp, 2)).Append('\t');
            }

            if (showLens)
                sb.Append(L(s.SeqLenA, 0)).Append('\t').Append(L(s.SeqLenB, 0)).Append('\t');
            if (showCov)
                sb.Append(CFormat.Fixed(Coverage(len1, s.SeqLenA), 2)).Append('\t').Append(CFormat.Fixed(Coverage(len2, s.SeqLenB), 2)).Append('\t');
            if (!isBrief && showDir)
                sb.Append(s.FrameA.ToString(CultureInfo.InvariantCulture)).Append('\t').Append(s.FrameB.ToString(CultureInfo.InvariantCulture)).Append('\t');
            sb.Append(s.IdA).Append('\t').Append(s.IdB);
            if (annotation)
                sb.Append('\t').Append(s.Annot);
            sb.Append('\n');
        }
    }
}
