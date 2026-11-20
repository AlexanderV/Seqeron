using System.Globalization;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-7 (B09 F172): evaluates one row of <c>emboss_seqopts_oracle.tsv</c> with the ports (what the EMBOSS program reported for the same input and qualifiers). A plain class so that the oracle test and a console
/// probe use the very same code.
/// </summary>
public static class EmbossSeqOptsProducer
{
    public sealed record Row(int Line, string Program, string Kind, string Params, string[] Seqs, string Quals1, string Quals2, string Expected)
    {
        public override string ToString() => $"line {Line} {Program} {Kind} [{Params}] {string.Join(" / ", Seqs.Where(s => s.Length > 0))} q1={Quals1} q2={Quals2}";
    }

    public static List<Row> ReadRows(TextReader reader)
    {
        var rows = new List<Row>();
        reader.ReadLine();
        int line = 1;
        string? text;
        while ((text = reader.ReadLine()) is not null)
        {
            line++;
            if (text.Length == 0)
                continue;
            var f = text.Split('\t');
            rows.Add(new Row(line, f[1], f[2], f[3], [f[4], f[5], f[6]], f[7], f[8], f[9]));
        }

        return rows;
    }

    public static EmbossSequenceQualifiers Quals(string text)
    {
        var f = text.Split('|').Select(x => int.Parse(x, CultureInfo.InvariantCulture)).ToArray();
        return new EmbossSequenceQualifiers { Begin = f[0], End = f[1], Reverse = f[2] == 1, Lower = f[3] == 1, Upper = f[4] == 1, Circular = f[5] == 1, Nucleotide = f[6] == 1, Protein = f[7] == 1 };
    }

    public static string Unescape(string s) => s.Replace("\\n", "\n", StringComparison.Ordinal).Replace("\\t", "\t", StringComparison.Ordinal).Replace("\\r", "\r", StringComparison.Ordinal).Replace("\\\\", "\\", StringComparison.Ordinal);

    /// <summary>What the program wrote for the row (unescaped; distmat's NaN sign is glibc's, the value is the same).</summary>
    public static string Expected(Row r)
    {
        if (r.Expected == "ERROR")
            return "ERROR";
        string text = Unescape(r.Expected);
        if (r.Program == "merger")
        {
            // merged sequence @@ the conflict lines of the report (position, base, position, base, used base)
            var parts = text.Split("@@", 2);
            var conflicts = System.Text.RegularExpressions.Regex.Matches(parts[1], @"#\s+(\d+)\s+'(.)'\s+(\d+)\s+'(.)'\s+'(.)'");
            return parts[0] + "@@" + string.Join("|", conflicts.Select(m => $"{m.Groups[1].Value}:{m.Groups[2].Value}:{m.Groups[3].Value}:{m.Groups[4].Value}:{m.Groups[5].Value}"));
        }

        if (r.Program == "megamerger")
        {
            var parts = text.Split("@@", 2);
            return parts[0].Length == 0 ? "" : parts[0] + "@@" + ReportLines(parts[1]);
        }
        return r.Program == "distmat" ? text.Replace("-nan", "nan", StringComparison.Ordinal) : text;
    }

    /// <summary>The port's answer in the form the oracle stores, or <c>ERROR</c> when the port refuses (the program refused the input).</summary>
    public static string Produce(Row r)
    {
        try
        {
            return ProduceCore(r);
        }
        catch (Exception e) when (e is ArgumentException or InvalidOperationException or NotSupportedException)
        {
            return "ERROR";
        }
    }

    private static string ProduceCore(Row r)
    {
        var q1 = Quals(r.Quals1);
        var q2 = Quals(r.Quals2);
        string a = r.Seqs[0], b = r.Seqs[1], c = r.Seqs[2];
        switch (r.Program)
        {
            case "needle":
            case "water":
            case "stretcher":
            case "matcher":
            case "supermatcher":
            case "needleall":
                return PairFile(r, q1, q2);
            case "merger":
                {
                    var m = SequenceAligner.MergerProgram(new EmbossSequence("a", a, "", q1), new EmbossSequence("b", b, "", q2));
                    return m.MergedSequence + "@@" + string.Join("|", m.Conflicts.Select(c => $"{c.PositionA}:{c.BaseA}:{c.PositionB}:{c.BaseB}:{c.Using}"));
                }
            case "megamerger":
                {
                    var m = GenomicAnalyzer.MegaMerge(new EmbossSequence("a", a, "", q1), new EmbossSequence("b", b, "", q2), int.Parse(r.Params, CultureInfo.InvariantCulture));
                    return m.MergedSequence.Length == 0 ? "" : m.MergedSequence + "@@" + ReportLines(m.Report);
                }
            case "fuzznuc":
                {
                    var p = r.Params.Split('|');
                    var hits = GenomicAnalyzer.FindPatternHits(a, q1, p[0], int.Parse(p[1], CultureInfo.InvariantCulture), p[2] == "1")
                        .Select(h => (Start: h.Position + 1, End: h.Position + h.Length, Strand: h.IsReverseStrand ? "-" : "+"))
                        .OrderBy(h => h.Start).ThenBy(h => h.End).ThenBy(h => h.Strand, StringComparer.Ordinal);
                    return Hits(hits.Select(h => $"{h.Start}:{h.End}:{h.Strand}"));
                }
            case "fuzzpro":
                {
                    var hits = GenomicAnalyzer.FindProteinPatternHits(a, q1, r.Params).OrderBy(h => h.Position).ThenBy(h => h.Length);
                    return Hits(hits.Select(h => $"{h.Position + 1}:{h.Position + h.Length}"));
                }
            case "fuzztran":
                {
                    var p = r.Params.Split('|');
                    var hits = GenomicAnalyzer.FindTranslatedPatternHits(a, q1, p[0], p[1])
                        .Select(h =>
                        {
                            // the gff report lists the smaller coordinate first, except when the two have opposite signs (a reverse-frame hit that straddles 0)
                            int s0 = h.EmbossStart!.Value, e0 = h.EmbossEnd!.Value;
                            bool swap = s0 > e0 && !(e0 < 0 && s0 >= 0);
                            return (Frame: Math.Abs(h.Frame), Start: swap ? e0 : s0, End: swap ? s0 : e0);
                        })
                        .OrderBy(h => h.Frame).ThenBy(h => h.Start).ThenBy(h => h.End);
                    return Hits(hits.Select(h => $"{h.Frame}:{h.Start}:{h.End}"));
                }
            case "etandem":
                {
                    var p = r.Params.Split(',').Select(x => int.Parse(x, CultureInfo.InvariantCulture)).ToArray();
                    var hits = GenomicAnalyzer.FindTandemRepeatsEtandem(a, q1, p[0], p[1], p[2], p[3] == 1, p[4] == 1);
                    return Hits2(hits.Select(h => string.Join(",", h.Score, h.Start, h.End, h.RepeatSize, h.RepeatCount, CFormat.Fixed(h.Identity, 1), h.Consensus)));
                }
            case "equicktandem":
                {
                    var p = r.Params.Split(',').Select(x => int.Parse(x, CultureInfo.InvariantCulture)).ToArray();
                    var hits = GenomicAnalyzer.FindTandemRepeatsEquicktandem(a, q1, p[0], p[1]);
                    return Hits2(hits.Select(h => string.Join(",", h.Score, h.Start, h.End, h.RepeatSize, h.RepeatCount)));
                }
            case "dotmatcher":
                {
                    var p = r.Params.Split(',');
                    var segs = DotMatcher.Compute(a, q1, b, q2, int.Parse(p[1], CultureInfo.InvariantCulture), int.Parse(p[2], CultureInfo.InvariantCulture), SubstitutionMatrix.Get(p[0]));
                    return Hits(segs.Select(s => $"{s.Start2}:{s.Start1}:{s.End2}:{s.End1}"));
                }
            case "seqmatchall":
                {
                    var seqs = r.Seqs.Where(x => x.Length > 0).ToList();
                    var matches = GenomicAnalyzer.SeqMatchAll(seqs, q1, int.Parse(r.Params, CultureInfo.InvariantCulture))
                        .OrderBy(m => m.Sequence1).ThenBy(m => m.Sequence2).ThenBy(m => m.Start1).ThenBy(m => m.End1).ThenBy(m => m.Start2).ThenBy(m => m.End2).ThenBy(m => m.Length);
                    return Hits(matches.Select(m => $"{m.Sequence1}:{m.Sequence2}:{m.Start1}..{m.End1}:{m.Start2}..{m.End2}:{m.Length}:{(m.Reversed ? "--" : "++")}"));
                }
            case "distmat":
                {
                    int method = int.Parse(r.Params, CultureInfo.InvariantCulture);
                    var rows = r.Seqs.Where(x => x.Length > 0).ToList();
                    var d = MultipleAlignmentStatistics.CalculateDistanceMatrix(rows, q1, new MsaDistanceOptions { Method = (MsaDistanceMethod)method }, (MsaDistanceMethod)(method <= 2 ? method : 0));
                    var cells = new List<string>();
                    for (int i = 0; i < rows.Count - 1; i++)
                        for (int j = i + 1; j < rows.Count; j++)
                            cells.Add(PrintedDistance(d[i, j]));
                    return string.Join(",", cells);
                }
            default:
                throw new InvalidOperationException(r.Program);
        }
    }

    /// <summary>distmat's print: D = D * (float)100, then <c>%.2f</c> (glibc: nan / inf).</summary>
    private static string PrintedDistance(double distance)
    {
        float v = (float)distance * 100f;
        if (float.IsNaN(v))
            return "nan";
        if (float.IsInfinity(v))
            return v > 0 ? "inf" : "-inf";
        return v.ToString("F2", CultureInfo.InvariantCulture);
    }

    private static string Hits2(IEnumerable<string> rows)
    {
        string s = string.Join("|", rows);
        return s.Length == 0 ? "-" : s;
    }

    // the report without comment lines and blank lines (the header is the program's run information)
    private static string ReportLines(string report) => string.Join("\n", report.Split('\n').Where(l => !l.StartsWith('#') && l.Trim().Length > 0).Select(l => l.TrimEnd()));

    private static string Hits(IEnumerable<string> hits)
    {
        string s = string.Join(" ", hits);
        return s.Length == 0 ? "-" : s;
    }

    private static string Body(string text)
    {
        string hashes = new('#', 40);
        var parts = text.Split(hashes + "\n");
        return parts.Length > 2 ? hashes + string.Join(hashes + "\n", parts.Skip(2)) : text;
    }

    private static string PairFile(Row r, EmbossSequenceQualifiers q1, EmbossSequenceQualifiers q2)
    {
        var first = new EmbossSequence("a", r.Seqs[0], "", q1);
        EmbossProgramResult res;
        double go = 0, ge = 0;
        bool integer = false, global = true;
        switch (r.Program)
        {
            case "needle":
            case "water":
                {
                    var g = r.Params.Split(',');
                    go = double.Parse(g[0], CultureInfo.InvariantCulture);
                    ge = double.Parse(g[1], CultureInfo.InvariantCulture);
                    var second = new EmbossSequence("b", r.Seqs[1], "", q2);
                    res = r.Program == "needle" ? SequenceAligner.NeedleProgram(first, second, null, go, ge) : SequenceAligner.WaterProgram(first, second, null, go, ge);
                    global = r.Program == "needle";
                    break;
                }
            case "stretcher":
                {
                    var g = r.Params.Split(',');
                    go = int.Parse(g[0], CultureInfo.InvariantCulture);
                    ge = int.Parse(g[1], CultureInfo.InvariantCulture);
                    res = SequenceAligner.StretcherProgram(first, new EmbossSequence("b", r.Seqs[1], "", q2), null, (int)go, (int)ge);
                    integer = true;
                    break;
                }
            case "matcher":
                global = false;
                integer = true;
                res = SequenceAligner.MatcherProgram(first, new EmbossSequence("b", r.Seqs[1], "", q2), null, null, 4, int.Parse(r.Params, CultureInfo.InvariantCulture));
                go = res.MatrixName == "EDNAFULL" ? 16 : 14;
                ge = 4;
                break;
            case "supermatcher":
            case "needleall":
                {
                    go = 10;
                    ge = 0.5;
                    global = r.Program == "needleall";
                    var targets = new[] { r.Seqs[1], r.Seqs[2] }.Where(x => x.Length > 0).Select((x, i) => new EmbossSequence("b" + (i + 1), x, "", q2)).ToList();
                    res = r.Program == "supermatcher"
                        ? SequenceAligner.Supermatcher([first], targets, null, 10, 0.5, int.Parse(r.Params, CultureInfo.InvariantCulture))
                        : SequenceAligner.NeedleAll([first], targets);
                    break;
                }
            default:
                throw new InvalidOperationException(r.Program);
        }

        var options = new EmbossAlignOptions
        {
            Program = r.Program,
            MatrixName = res.MatrixName,
            Matrix = SubstitutionMatrix.Get(res.MatrixName),
            Global = global,
            IntegerScores = integer,
            StatisticsSetByProgram = integer,
            GapPenalty = go,
            ExtendPenalty = ge,
        };
        string text = EmbossAlignmentWriter.Format(EmbossAlignFormat.Pair, res.Alignments, options);
        return Body(text);
    }
}
