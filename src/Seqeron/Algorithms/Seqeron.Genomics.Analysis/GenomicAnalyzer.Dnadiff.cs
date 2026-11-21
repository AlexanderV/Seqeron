using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>The files <c>dnadiff</c> writes (<c>out.report</c>, <c>out.delta</c>, <c>.1delta</c>, <c>.mdelta</c>, <c>.1coords</c>, <c>.mcoords</c>, <c>.snps</c>, <c>.rdiff</c>, <c>.qdiff</c>, <c>.unref</c>, <c>.unqry</c>).</summary>
/// <param name="Report">Summary of alignments, differences and SNPs.</param>
/// <param name="Delta">Standard nucmer alignment output (<c>nucmer --maxmatch</c>).</param>
/// <param name="Delta1">1-to-1 alignment (<c>delta-filter -1</c>).</param>
/// <param name="DeltaM">M-to-M alignment (<c>delta-filter -m</c>).</param>
/// <param name="Coords1">1-to-1 coordinates (<c>show-coords -rclTH</c>).</param>
/// <param name="CoordsM">M-to-M coordinates (<c>show-coords -rclTH</c>).</param>
/// <param name="Snps">SNPs (<c>show-snps -rlTHC</c> of the 1-to-1 alignment).</param>
/// <param name="ReferenceDiff">Classified reference breakpoints (<c>show-diff -rH</c>).</param>
/// <param name="QueryDiff">Classified query breakpoints (<c>show-diff -qH</c>).</param>
/// <param name="UnalignedReference">Unaligned reference IDs and lengths (null when every reference sequence aligned).</param>
/// <param name="UnalignedQuery">Unaligned query IDs and lengths (null when every query sequence aligned).</param>
public sealed record DnadiffResult(string Report, string Delta, string Delta1, string DeltaM, string Coords1, string CoordsM, string Snps, string ReferenceDiff, string QueryDiff, string? UnalignedReference, string? UnalignedQuery);

public static partial class GenomicAnalyzer
{
    #region dnadiff

    private const int DnadiffSnpBuffer = 20;

    /// <summary>
    /// MUMmer 3.23 <c>dnadiff</c> (a perl script over <c>nucmer --maxmatch</c>, <c>delta-filter -1 / -m</c>, <c>show-coords -rclTH</c>, <c>show-snps -rlTHC</c>, <c>show-diff -rH / -qH</c>): the comparative analysis of two sequence sets
    /// and its report (sequences, bases, 1-to-1 and M-to-M alignments, breakpoints, relocations, translocations, inversions, insertions, tandem insertions, SNPs and indels with the 12 substitution classes, all
    /// also for the "good" ones whose buffer is at least 20).
    /// </summary>
    /// <remarks>
    /// The order of the per-substitution lines inside the SNP / indel blocks of the report and the order of the <c>.unref</c> / <c>.unqry</c> lines is the order of perl's hash iteration, which differs from run to run
    /// (hash randomisation; measured: repeated runs of the installed program differ); here they are in a fixed order (. A C G T, then any other symbol). Everything else is byte for byte the program's.
    /// </remarks>
    /// <param name="referenceFasta">Reference FASTA text.</param>
    /// <param name="queryFasta">Query FASTA text.</param>
    /// <param name="referenceName">Reference file name for the delta header (dnadiff passes the absolute path).</param>
    /// <param name="queryName">Query file name for the delta header.</param>
    /// <param name="deltaText">A precomputed delta file (<c>dnadiff -d</c>); null = run <c>nucmer --maxmatch</c>.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="FormatException">A program aborts on the input (the perl script's "Whitespace found in FastA", "Unrecognized format").</exception>
    public static DnadiffResult Dnadiff(string referenceFasta, string queryFasta, string referenceName = "reference.fa", string queryName = "query.fa", string? deltaText = null, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(referenceFasta);
        ArgumentNullException.ThrowIfNull(queryFasta);
        string delta = deltaText ?? Nucmer(referenceFasta, queryFasta, referenceName, queryName, new NucmerOptions { MatchMode = CommonRegionMode.MaximalExactMatches }, cancellationToken).Delta!;
        string delta1 = MummerDeltaFilter.Run(delta, new MummerDeltaFilterOptions { OneToOne = true });
        string deltaM = MummerDeltaFilter.Run(delta, new MummerDeltaFilterOptions { ManyToMany = true });
        var coordsOptions = new MummerShowCoordsOptions { SortByReference = true, ShowCoverage = true, ShowSeqLengths = true, Tabular = true, PrintHeader = false };
        string coords1 = MummerShowCoords.Run(delta1, coordsOptions).Text;
        string coordsM = MummerShowCoords.Run(deltaM, coordsOptions).Text;
        string snps = MummerShowSnps.Run(delta1, referenceFasta, queryFasta, new MummerShowSnpsOptions { SortReference = true, ShowLength = true, Tabular = true, PrintHeader = false, ShowConflict = false }).Text;
        string rdiff = MummerShowDiff.Run(deltaM, new MummerShowDiffOptions { ReferenceDiff = true, PrintHeader = false }).Text;
        string qdiff = MummerShowDiff.Run(deltaM, new MummerShowDiffOptions { QueryDiff = true, PrintHeader = false }).Text;
        return DnadiffReport(delta, coords1, coordsM, snps, rdiff, qdiff, referenceFasta, queryFasta, delta1, deltaM);
    }

    private static Dictionary<string, long> FastaSizes(string fasta, List<string> order)
    {
        var sizes = new Dictionary<string, long>(StringComparer.Ordinal);
        string? tag = null;
        long len = 0;
        foreach (string raw in fasta.Split('\n'))
        {
            string line = raw;
            if (line.StartsWith('>'))
            {
                if (tag != null)
                {
                    sizes[tag] = len;
                    if (!order.Contains(tag))
                        order.Add(tag);
                }

                var m = System.Text.RegularExpressions.Regex.Match(line, @"^>(\S+)");
                tag = m.Success ? m.Groups[1].Value : null;
                len = 0;
            }
            else
            {
                if (line.Any(char.IsWhiteSpace))
                    throw new FormatException("ERROR: Whitespace found in FastA, aborting.");
                len += line.Length;
            }
        }

        if (tag != null)
        {
            sizes[tag] = len;
            if (!order.Contains(tag))
                order.Add(tag);
        }

        return sizes;
    }

    private static DnadiffResult DnadiffReport(string delta, string coords1, string coordsM, string snps, string rdiff, string qdiff, string refFasta, string qryFasta, string delta1, string deltaM)
    {
        var ci = CultureInfo.InvariantCulture;
        var lines = delta.Split('\n');
        string header = lines[0] + "\n" + lines[1] + "\n" + "\n";
        var refOrder = new List<string>();
        var qryOrder = new List<string>();
        var refs = FastaSizes(refFasta, refOrder);
        var qrys = FastaSizes(qryFasta, qryOrder);

        long rnSeqs = 0, qnSeqs = 0, rnBases = 0, qnBases = 0, rnASeqs = 0, qnASeqs = 0, rnABases = 0, qnABases = 0, rnBrk = 0, qnBrk = 0;
        foreach (long v in refs.Values)
        {
            rnSeqs++;
            rnBases += v;
        }

        foreach (long v in qrys.Values)
        {
            qnSeqs++;
            qnBases += v;
        }

        long rqnAlignsM = 0, rSumLenM = 0, qSumLenM = 0;
        double rqSumIdyM = 0, rqSumLenM = 0;
        foreach (string line in coordsM.Split('\n'))
        {
            if (line.Length == 0)
                continue;
            var a = line.Split('\t');
            if (a.Length != 13)
                throw new FormatException("ERROR: Unrecognized format .mcoords, aborting.");
            rqnAlignsM++;
            rSumLenM += long.Parse(a[4], ci);
            qSumLenM += long.Parse(a[5], ci);
            rqSumIdyM += (double.Parse(a[6], ci) / 100.0) * (double.Parse(a[4], ci) + double.Parse(a[5], ci));
            rqSumLenM += double.Parse(a[4], ci) + double.Parse(a[5], ci);
            if (refs[a[11]] > 0)
            {
                rnASeqs++;
                rnABases += refs[a[11]];
                refs[a[11]] *= -1;
            }

            if (qrys[a[12]] > 0)
            {
                qnASeqs++;
                qnABases += qrys[a[12]];
                qrys[a[12]] *= -1;
            }

            long lo, hi, s0 = long.Parse(a[0], ci), s1 = long.Parse(a[1], ci), s2 = long.Parse(a[2], ci), s3 = long.Parse(a[3], ci);
            (lo, hi) = s0 < s1 ? (s0, s1) : (s1, s0);
            if (lo != 1)
                rnBrk++;
            if (hi != long.Parse(a[7], ci))
                rnBrk++;
            (lo, hi) = s2 < s3 ? (s2, s3) : (s3, s2);
            if (lo != 1)
                qnBrk++;
            if (hi != long.Parse(a[8], ci))
                qnBrk++;
        }

        long rqnAligns1 = 0, rSumLen1 = 0, qSumLen1 = 0;
        double rqSumIdy1 = 0, rqSumLen1 = 0;
        foreach (string line in coords1.Split('\n'))
        {
            if (line.Length == 0)
                continue;
            var a = line.Split('\t');
            if (a.Length != 13)
                throw new FormatException("ERROR: Unrecognized format .1coords, aborting.");
            rqnAligns1++;
            rSumLen1 += long.Parse(a[4], ci);
            qSumLen1 += long.Parse(a[5], ci);
            rqSumIdy1 += (double.Parse(a[6], ci) / 100.0) * (double.Parse(a[4], ci) + double.Parse(a[5], ci));
            rqSumLen1 += double.Parse(a[4], ci) + double.Parse(a[5], ci);
        }

        long rnIns = 0, qnIns = 0, rSumIns = 0, qSumIns = 0, rnTIns = 0, qnTIns = 0, rSumTIns = 0, qSumTIns = 0, rnInv = 0, qnInv = 0, rnRel = 0, qnRel = 0, rnTrn = 0, qnTrn = 0;
        void Diffs(string text, bool reference)
        {
            foreach (string line in text.Split('\n'))
            {
                if (line.Length == 0)
                    continue;
                var a = line.Split('\t');
                if (a.Length < 5)
                    throw new FormatException("ERROR: Unrecognized format diff, aborting.");
                long gap = long.Parse(a[4], ci);
                long ins = gap;
                if (a[1] == "GAP")
                {
                    if (a.Length != 7)
                        throw new FormatException("ERROR: Unrecognized format diff, aborting.");
                    long a5 = long.Parse(a[5], ci), a6 = long.Parse(a[6], ci);
                    if (a6 > gap)
                        ins = a6;
                    if (gap <= 0 && a5 <= 0 && a6 > 0)
                    {
                        if (reference)
                        {
                            rnTIns++;
                            rSumTIns += a6;
                        }
                        else
                        {
                            qnTIns++;
                            qSumTIns += a6;
                        }
                    }
                }

                if (a[1] != "DUP" && gap > 0)
                {
                    if (reference)
                        rnABases -= gap;
                    else
                        qnABases -= gap;
                }

                if (ins > 0)
                {
                    if (reference)
                    {
                        rnIns++;
                        rSumIns += ins;
                    }
                    else
                    {
                        qnIns++;
                        qSumIns += ins;
                    }
                }

                if (a[1] == "INV")
                    if (reference) rnInv++; else qnInv++;
                if (a[1] == "JMP")
                    if (reference) rnRel++; else qnRel++;
                if (a[1] == "SEQ")
                    if (reference) rnTrn++; else qnTrn++;
            }
        }

        Diffs(rdiff, true);
        Diffs(qdiff, false);

        var rqSnps = NewSubstitutionTable();
        var rqGSnps = NewSubstitutionTable();
        long rqnSnps = 0, rqnIndels = 0, rqnGSnps = 0, rqnGIndels = 0;
        foreach (string line in snps.Split('\n'))
        {
            if (line.Length == 0)
                continue;
            var a = line.Split('\t');
            if (a.Length != 12)
                throw new FormatException("ERROR: Unrecognized format .snps, aborting");
            string r = a[1].ToUpperInvariant(), q = a[2].ToUpperInvariant();
            Bump(rqSnps, r, q);
            if (r == "." || q == ".")
                rqnIndels++;
            else
                rqnSnps++;
            if (long.Parse(a[4], ci) >= DnadiffSnpBuffer)
            {
                Bump(rqGSnps, r, q);
                if (r == "." || q == ".")
                    rqnGIndels++;
                else
                    rqnGSnps++;
            }
        }

        var sb = new StringBuilder();
        string F(string label, string x, string y) => label.PadRight(15) + " " + x.PadLeft(20) + " " + y.PadLeft(20) + "\n";
        string D(string label, long x, long y) => F(label, x.ToString(ci), y.ToString(ci));
        string Fx(string label, double x, double y) => F(label, CFormat.Fixed(x, 2), CFormat.Fixed(y, 2));
        string Pc(long n, double pct) => n.ToString(ci).PadLeft(10) + "(" + CFormat.Fixed(pct, 2) + "%)";
        double Div(double a, double b) => a / b;

        sb.Append(header);
        sb.Append(F(string.Empty, "[REF]", "[QRY]"));
        sb.Append("[Sequences]\n");
        sb.Append(D("TotalSeqs", rnSeqs, qnSeqs));
        sb.Append(F("AlignedSeqs", Pc(rnASeqs, rnSeqs != 0 ? Div(rnASeqs, rnSeqs) * 100.0 : 0), Pc(qnASeqs, rnSeqs != 0 ? Div(qnASeqs, qnSeqs) * 100.0 : 0)));
        sb.Append(F("UnalignedSeqs", Pc(rnSeqs - rnASeqs, rnSeqs != 0 ? Div(rnSeqs - rnASeqs, rnSeqs) * 100.0 : 0), Pc(qnSeqs - qnASeqs, qnSeqs != 0 ? Div(qnSeqs - qnASeqs, qnSeqs) * 100.0 : 0)));
        sb.Append("\n[Bases]\n");
        sb.Append(D("TotalBases", rnBases, qnBases));
        sb.Append(F("AlignedBases", Pc(rnABases, rnBases != 0 ? Div(rnABases, rnBases) * 100.0 : 0), Pc(qnABases, qnBases != 0 ? Div(qnABases, qnBases) * 100.0 : 0)));
        sb.Append(F("UnalignedBases", Pc(rnBases - rnABases, rnBases != 0 ? Div(rnBases - rnABases, rnBases) * 100.0 : 0), Pc(qnBases - qnABases, qnBases != 0 ? Div(qnBases - qnABases, qnBases) * 100.0 : 0)));
        sb.Append("\n[Alignments]\n");
        sb.Append(D("1-to-1", rqnAligns1, rqnAligns1));
        sb.Append(D("TotalLength", rSumLen1, qSumLen1));
        sb.Append(Fx("AvgLength", rqnAligns1 != 0 ? Div(rSumLen1, rqnAligns1) : 0, rqnAligns1 != 0 ? Div(qSumLen1, rqnAligns1) : 0));
        double idy1 = rqSumLen1 != 0 ? Div(rqSumIdy1, rqSumLen1) * 100.0 : 0;
        sb.Append(Fx("AvgIdentity", idy1, idy1));
        sb.Append('\n');
        sb.Append(D("M-to-M", rqnAlignsM, rqnAlignsM));
        sb.Append(D("TotalLength", rSumLenM, qSumLenM));
        sb.Append(Fx("AvgLength", rqnAlignsM != 0 ? Div(rSumLenM, rqnAlignsM) : 0, rqnAlignsM != 0 ? Div(qSumLenM, rqnAlignsM) : 0));
        double idyM = rqSumLenM != 0 ? Div(rqSumIdyM, rqSumLenM) * 100.0 : 0;
        sb.Append(Fx("AvgIdentity", idyM, idyM));
        sb.Append("\n[Feature Estimates]\n");
        sb.Append(D("Breakpoints", rnBrk, qnBrk));
        sb.Append(D("Relocations", rnRel, qnRel));
        sb.Append(D("Translocations", rnTrn, qnTrn));
        sb.Append(D("Inversions", rnInv, qnInv));
        sb.Append('\n');
        sb.Append(D("Insertions", rnIns, qnIns));
        sb.Append(D("InsertionSum", rSumIns, qSumIns));
        sb.Append(Fx("InsertionAvg", rnIns != 0 ? Div(rSumIns, rnIns) : 0, qnIns != 0 ? Div(qSumIns, qnIns) : 0));
        sb.Append('\n');
        sb.Append(D("TandemIns", rnTIns, qnTIns));
        sb.Append(D("TandemInsSum", rSumTIns, qSumTIns));
        sb.Append(Fx("TandemInsAvg", rnTIns != 0 ? Div(rSumTIns, rnTIns) : 0, qnTIns != 0 ? Div(qSumTIns, qnTIns) : 0));
        sb.Append("\n[SNPs]\n");
        sb.Append(D("TotalSNPs", rqnSnps, rqnSnps));
        SubstitutionLines(sb, rqSnps, rqnSnps, (r, q) => r != "." && q != ".", ci, F, Pc);
        sb.Append('\n');
        sb.Append(D("TotalGSNPs", rqnGSnps, rqnGSnps));
        SubstitutionLines(sb, rqGSnps, rqnGSnps, (r, q) => r != "." && q != ".", ci, F, Pc);
        sb.Append('\n');
        sb.Append(D("TotalIndels", rqnIndels, rqnIndels));
        SubstitutionLines(sb, rqSnps, rqnIndels, (r, q) => q == ".", ci, F, Pc);
        SubstitutionLines(sb, rqSnps, rqnIndels, (r, q) => r == ".", ci, F, Pc);
        sb.Append('\n');
        sb.Append(D("TotalGIndels", rqnGIndels, rqnGIndels));
        SubstitutionLines(sb, rqGSnps, rqnGIndels, (r, q) => q == ".", ci, F, Pc);
        SubstitutionLines(sb, rqGSnps, rqnGIndels, (r, q) => r == ".", ci, F, Pc);

        string? unref = null, unqry = null;
        if (rnSeqs != rnASeqs)
            unref = string.Concat(refOrder.Where(k => refs[k] >= 0).Select(k => $"{k}\tUNI\t1\t{refs[k]}\t{refs[k]}\n"));
        if (qnSeqs != qnASeqs)
            unqry = string.Concat(qryOrder.Where(k => qrys[k] >= 0).Select(k => $"{k}\tUNI\t1\t{qrys[k]}\t{qrys[k]}\n"));
        return new DnadiffResult(sb.ToString(), delta, delta1, deltaM, coords1, coordsM, snps, rdiff, qdiff, unref, unqry);
    }

    private static SortedDictionary<string, SortedDictionary<string, long>> NewSubstitutionTable()
    {
        var order = new[] { ".", "A", "C", "G", "T" };
        var t = new SortedDictionary<string, SortedDictionary<string, long>>(Comparer<string>.Create(CompareSymbols));
        foreach (string r in order)
        {
            t[r] = new SortedDictionary<string, long>(Comparer<string>.Create(CompareSymbols));
            foreach (string q in order)
                if (q != r)
                    t[r][q] = 0;
        }

        return t;
    }

    // . A C G T first (the keys of the program's hash), then any other symbol in ordinal order
    private static int CompareSymbols(string a, string b)
    {
        int Rank(string s) => s switch { "." => 0, "A" => 1, "C" => 2, "G" => 3, "T" => 4, _ => 5 };
        int c = Rank(a).CompareTo(Rank(b));
        return c != 0 ? c : string.CompareOrdinal(a, b);
    }

    private static void Bump(SortedDictionary<string, SortedDictionary<string, long>> t, string r, string q)
    {
        if (!t.TryGetValue(r, out var row))
            t[r] = row = new SortedDictionary<string, long>(Comparer<string>.Create(CompareSymbols));
        row[q] = row.GetValueOrDefault(q) + 1;
        if (!t.TryGetValue(q, out var rowQ))
            t[q] = rowQ = new SortedDictionary<string, long>(Comparer<string>.Create(CompareSymbols));
        if (!rowQ.ContainsKey(r))
            rowQ[r] = 0;
    }

    private static void SubstitutionLines(StringBuilder sb, SortedDictionary<string, SortedDictionary<string, long>> t, long total, Func<string, string, bool> select, CultureInfo ci,
        Func<string, string, string, string> line, Func<long, double, string> pc)
    {
        foreach (var (r, row) in t)
        {
            foreach (var (q, n) in row)
            {
                if (!select(r, q))
                    continue;
                long other = t[q].GetValueOrDefault(r);
                sb.Append(line(r + q, pc(n, total != 0 ? (double)n / total * 100.0 : 0), pc(other, total != 0 ? (double)other / total * 100.0 : 0)));
            }
        }
    }

    #endregion
}
