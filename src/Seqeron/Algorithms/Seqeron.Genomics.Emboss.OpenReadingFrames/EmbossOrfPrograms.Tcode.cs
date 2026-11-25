// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.OpenReadingFrames, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127, S1905
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.OpenReadingFrames;

using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Emboss.SequenceFormats.EmbossSequencePrograms;

/// <summary><c>tcode</c> of EMBOSS 6.6.0 (<c>emboss/tcode.c</c>): the Fickett TESTCODE statistic of a nucleotide sequence in a sliding window, as a report (B09 WP63, F282).</summary>
public static partial class EmbossOrfPrograms
{
    /// <summary>The data file <c>Etcode.dat</c> of EMBOSS 6.6.0 (the default of <c>-datafile</c>).</summary>
    internal const string EtcodeDat = """
# Fickett TESTCODE data
# NAR 10(17) 5303-5318
#
# Position parameter values
1.9
1.8
1.7
1.6
1.5
1.4
1.3
1.2
1.1
0.0
#
#
# Content parameter values
0.33
0.31
0.29
0.27
0.25
0.23
0.21
0.19
0.17
0.00
#
#
# Position probabilities for A,C,G,T respectively
0.94 0.80 0.90 0.97
0.68 0.70 0.88 0.97
0.84 0.70 0.74 0.91
0.93 0.81 0.64 0.68
0.58 0.66 0.53 0.69
0.68 0.48 0.48 0.44
0.45 0.51 0.27 0.54
0.34 0.33 0.16 0.20
0.20 0.30 0.08 0.09
0.22 0.23 0.08 0.09
#
#
# Content probabilities for A,C,G,T respectively
0.28 0.82 0.40 0.28
0.49 0.64 0.54 0.24
0.44 0.51 0.47 0.39
0.55 0.64 0.64 0.40
0.62 0.59 0.64 0.55
0.49 0.59 0.73 0.75
0.67 0.43 0.41 0.56
0.65 0.44 0.41 0.69
0.81 0.39 0.33 0.51
0.21 0.31 0.29 0.58
#
#
# Weights for position
0.26
0.18
0.31
0.33
#
#
# Weights for content
0.11
0.12
0.15
0.14
""";

    private sealed class TestcodeTable
    {
        public readonly List<float> Positions = [], Content = [], PprobA = [], PprobC = [], PprobG = [], PprobT = [], CprobA = [], CprobC = [], CprobG = [], CprobT = [], Pweights = [], Cweights = [];
    }

    // ajFloatGet of the array: an index beyond the data reads memory the library never wrote (the program has no defined result)
    private static float TcodeGet(List<float> a, int i)
        => i >= 0 && i < a.Count ? a[i] : throw new InvalidOperationException("tcode reads an element beyond the end of a data array of the testcode file (ajFloatGet: bad float array index): the reference program has no defined result.");

    // ajFmtScanS(line, "%f%f%f%f", ...): the values that are scanned are replaced, a value that does not scan keeps the one before
    private static void TcodeScan(string line, float[] vals, int count)
    {
        int pos = 0;
        for (int k = 0; k < count; k++)
        {
            while (pos < line.Length && IsCSpace(line[pos]))
                pos++;
            int start = pos;
            if (pos < line.Length && (line[pos] == '+' || line[pos] == '-'))
                pos++;
            int digits = 0;
            while (pos < line.Length && char.IsAsciiDigit(line[pos]))
            {
                pos++;
                digits++;
            }

            if (pos < line.Length && line[pos] == '.')
            {
                pos++;
                while (pos < line.Length && char.IsAsciiDigit(line[pos]))
                {
                    pos++;
                    digits++;
                }
            }

            if (digits == 0)
                return;
            if (pos < line.Length && (line[pos] == 'e' || line[pos] == 'E'))
            {
                int save = pos++;
                if (pos < line.Length && (line[pos] == '+' || line[pos] == '-'))
                    pos++;
                int ed = 0;
                while (pos < line.Length && char.IsAsciiDigit(line[pos]))
                {
                    pos++;
                    ed++;
                }

                if (ed == 0)
                    pos = save;
            }

            vals[k] = float.Parse(line.AsSpan(start, pos - start), NumberStyles.Float, CultureInfo.InvariantCulture);
        }
    }

    private static bool IsCSpace(char c) => c is ' ' or '\t' or '\n' or '\v' or '\f' or '\r';

    // tcode_readdata: the sections of the data file are read line by line (ajReadlineTrim: the line without its white space at the end; a line that is blank is not a comment)
    private static bool TcodeReadData(TestcodeTable t, string text)
    {
        var lines = new List<string>();
        int p = 0;
        while (p < text.Length)
        {
            int nl = text.IndexOf('\n', p);
            string line = nl < 0 ? text[p..] : text[p..nl];
            p = nl < 0 ? text.Length : nl + 1;
            lines.Add(line.TrimEnd(' ', '\t', '\r', '\n', '\v', '\f'));
        }

        int li = 0;
        bool Read(out string line)
        {
            if (li < lines.Count)
            {
                line = lines[li++];
                return true;
            }

            line = "";
            return false;
        }

        bool ok = true;
        var v = new float[4];
        v[0] = float.MaxValue;                                // BIGVAL
        while (v[0] != 0 && ok)
        {
            ok = Read(out string line);
            if (line.Length > 0 && (line[0] == '#' || line[0] == '\n'))
                continue;
            TcodeScan(line, v, 1);
            t.Positions.Add(v[0]);
        }

        ok = true;
        v[0] = float.MaxValue;
        while (v[0] != 0 && ok)
        {
            ok = Read(out string line);
            if (line.Length > 0 && (line[0] == '#' || line[0] == '\n'))
                continue;
            TcodeScan(line, v, 1);
            t.Content.Add(v[0]);
        }

        for (int pass = 0; pass < 2; pass++)
        {
            int n = pass == 0 ? t.Positions.Count : t.Content.Count;
            int i = 0;
            ok = true;
            while (ok && i != n)
            {
                ok = Read(out string line);
                if (line.Length > 0 && (line[0] == '#' || line[0] == '\n'))
                    continue;
                TcodeScan(line, v, 4);
                (pass == 0 ? t.PprobA : t.CprobA).Add(v[0]);
                (pass == 0 ? t.PprobC : t.CprobC).Add(v[1]);
                (pass == 0 ? t.PprobG : t.CprobG).Add(v[2]);
                (pass == 0 ? t.PprobT : t.CprobT).Add(v[3]);
                i++;
            }
        }

        for (int pass = 0; pass < 2; pass++)
        {
            int i = 0;
            ok = true;
            while (ok && i != 4)
            {
                ok = Read(out string line);
                if (line.Length > 0 && (line[0] == '#' || line[0] == '\n'))
                    continue;
                TcodeScan(line, v, 1);
                (pass == 0 ? t.Pweights : t.Cweights).Add(v[0]);
                i++;
            }
        }

        return ok;
    }

    private static int TcodeIndex(List<float> array, float value)
    {
        if (float.IsNaN(value))
            throw new InvalidOperationException("tcode_index never ends on a window without A, C, G or T (the composition is not a number): the reference program loops forever.");
        int i = 0;
        while (true)
        {
            float thisval = TcodeGet(array, i);
            if (value >= thisval)
                return i;
            i++;
        }
    }

    // tcode_slide
    private static float TcodeSlide(string p, int window, TestcodeTable table, int pos)
    {
        var scores = new int[3, 4];
        for (int i = 0; i < window; ++i)
        {
            switch (p[i + pos])
            {
                case 'A':
                    ++scores[i % 3, 0];
                    break;
                case 'C':
                    ++scores[i % 3, 1];
                    break;
                case 'G':
                    ++scores[i % 3, 2];
                    break;
                case 'T':
                    ++scores[i % 3, 3];
                    break;
            }
        }

        int asum = scores[0, 0] + scores[1, 0] + scores[2, 0];
        int csum = scores[0, 1] + scores[1, 1] + scores[2, 1];
        int gsum = scores[0, 2] + scores[1, 2] + scores[2, 2];
        int tsum = scores[0, 3] + scores[1, 3] + scores[2, 3];
        int sum = asum + csum + gsum + tsum;
        float acomp = (float)asum / (float)sum, ccomp = (float)csum / (float)sum, gcomp = (float)gsum / (float)sum, tcomp = (float)tsum / (float)sum;
        float Pos(int col)
        {
            int vmax = Math.Max(Math.Max(scores[0, col], scores[1, col]), scores[2, col]);
            int vmin = Math.Min(Math.Min(scores[0, col], scores[1, col]), scores[2, col]);
            return (float)((float)vmax / (float)((float)vmin + 1.0));
        }

        float apos = Pos(0), cpos = Pos(1), gpos = Pos(2), tpos = Pos(3);
        float p1 = (float)(TcodeGet(table.PprobA, TcodeIndex(table.Positions, apos)) * TcodeGet(table.Pweights, 0));
        float p2 = (float)(TcodeGet(table.PprobC, TcodeIndex(table.Positions, cpos)) * TcodeGet(table.Pweights, 1));
        float p3 = (float)(TcodeGet(table.PprobG, TcodeIndex(table.Positions, gpos)) * TcodeGet(table.Pweights, 2));
        float p4 = (float)(TcodeGet(table.PprobT, TcodeIndex(table.Positions, tpos)) * TcodeGet(table.Pweights, 3));
        float p5 = (float)(TcodeGet(table.CprobA, TcodeIndex(table.Content, acomp)) * TcodeGet(table.Cweights, 0));
        float p6 = (float)(TcodeGet(table.CprobC, TcodeIndex(table.Content, ccomp)) * TcodeGet(table.Cweights, 1));
        float p7 = (float)(TcodeGet(table.CprobG, TcodeIndex(table.Content, gcomp)) * TcodeGet(table.Cweights, 2));
        float p8 = (float)(TcodeGet(table.CprobT, TcodeIndex(table.Content, tcomp)) * TcodeGet(table.Cweights, 3));
        float r = (float)(p1 + p2);
        r = (float)(r + p3);
        r = (float)(r + p4);
        r = (float)(r + p5);
        r = (float)(r + p6);
        r = (float)(r + p7);
        return (float)(r + p8);
    }

    /// <summary>
    /// <c>tcode</c>: for every sequence the TESTCODE statistic of the windows of <c>-window</c> bases every <c>-step</c> bases (<c>tcode_slide</c> with the probability tables of <c>-datafile</c>), reported as the table of
    /// <c>-outfile</c> (<c>-rformat</c> table by default; an estimation Non-coding below 0.74, Coding from 0.95, else No opinion); <c>-plot</c> draws the values (<c>-graph data</c> / <c>none</c>).
    /// </summary>
    internal static EmbossRunResult Tcode(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("tcode");
        var all = EmbossProgramBase.Open("tcode", a, files, ctx, "dna");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        var first = all.Peek();
        string dataText;
        int window, step;
        bool plot;
        EmbossReportArgs.Opened? report = null;
        EmbossGraphXy? graph;
        try
        {
            string? dataName = a.Raw("datafile");
            if (dataName is null || (dataName == "Etcode.dat" && !files.ContainsKey(dataName)))
                dataText = EtcodeDat;
            else if (dataName.Length == 0)
            {
                ctx.Error("Input file is required");
                throw a.Die("datafile");
            }
            else if (files.TryGetValue(dataName, out string? dt))
                dataText = dt;
            else
            {
                ctx.Error($"Unable to open data file '{dataName}' for input");
                throw a.Die("datafile");
            }

            int len0 = first.Seq.Length;
            window = a.IntClamped("window", 200, 200, len0, ctx);
            step = a.IntClamped("step", 3, 1, len0 - window, ctx);
            plot = a.Bool("plot", false);
            string? given = a.Raw("outfile");
            _ = EmbossReportArgs.Qualifiers(a, "outfile");                      // the values are read whether the report is written or not
            if (plot && given is { Length: 0 })
                given = null;                                                   // nullok + nulldefault: an empty value makes the default name
            if (!(plot && given is null && a.Raw("outfile") is null))
            {
                report = EmbossReportArgs.Open(a, ctx, "outfile", "tcode", "table", [new("str", "Estimation", "Estimation")], 3, "Fickett TESTCODE statistic", given,
                    EmbossOutPath.DefaultBaseName(first.Name, a.Raw("sequence") ?? ""), EmbossRunEnvironment.Current.RunDate ?? DateTime.Now, files, multi: true);
                if (report is null)
                    return EmbossProgramBase.Died(ctx, "outfile");
            }

            graph = EmbossGraphArgs.Open(a, ctx, "tcode", "graph", "Fickett testcode plot", 1, EmbossOutPath.DefaultBaseName(first.Name, a.Raw("sequence") ?? ""), files);
            if (graph is null)
            {
                var died = EmbossProgramBase.Died(ctx, "graph");
                return report is null || report.File == "stdout" ? died : died with { Files = new Dictionary<string, string> { [report.File] = "" } };
            }
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var table = new TestcodeTable();
        if (!TcodeReadData(table, dataText))
            return EmbossProgramBase.Fatal(ctx, "\n   EMBOSS An error in tcode.c at line 164:\nIncorrect format Testcode data file");
        var tables = new List<EmbossReport>();
        long count = 0, total = 0;
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            count++;
            EmbossSeqRange.Trim(seq);
            total += seq.Seq.Length;
            string seqstr = seq.SeqText.ToUpperInvariant();
            AssignUpper(seq, seqstr);                                          // ajSeqFmtUpper: the sequence object itself is in upper case (the report prints it)
            int len = seqstr.Length;
            if (len - window + 1 < 0)
                throw new InvalidOperationException("tcode puts the array element len - window + 1 below 0 (ajFloatPut of an unsigned index: 'Memory allocation failed in ajMemCallocZero'): the reference program has no defined result.");
            int pos = 0;
            var from = new List<int>();
            var to = new List<int>();
            var codes = new List<float>();
            while (pos + window < len)
            {
                if (step == 0)
                    throw new InvalidOperationException("tcode with a step of 0 never leaves the first window (the memory of its arrays grows until the program is killed): the reference program has no defined result.");
                float tc = TcodeSlide(seqstr, window, table, pos);
                from.Add(pos + 1 + seq.Offset);
                to.Add(pos + window + seq.Offset);
                codes.Add(tc);
                pos += step;
            }

            if (report is not null)
                tables.Add(TcodeReport(seq, from, to, codes, count == 1));
            if (plot)
                TcodePlot(graph, seq, from, to, codes, ctx);
        }

        string text = report is null ? "" : EmbossReportArgs.Write(report, tables, ctx, count, total);
        string created = graph.Close();                                         // ajGraphicsClose: "Created tcode1.dat" on the standard output
        if (report is null)
            return new EmbossRunResult(created, EmbossRunResult.Join(ctx.Messages), 0, new Dictionary<string, string>(graph.Files));
        var result = EmbossReportArgs.Result(ctx, report.File, text, graph.Files);
        return result with { Stdout = result.Stdout + created };
    }

    // the plot of one sequence: the mean position of every window against its statistic, with the lines at 0.74 (red) and 0.95 (green)
    private static void TcodePlot(EmbossGraphXy graph, EmbossSeq seq, List<int> from, List<int> to, List<float> codes, EmbossFeatureReadContext ctx)
    {
        int n = codes.Count;
        if (n == 0)
            throw new InvalidOperationException("tcode plots a sequence without a window: x[-1] is read before the array (valgrind: invalid read): the reference program has no defined result.");
        var g = new EmbossGraphData { Type = "2D plot", Dataname = seq.Name, X = new float[n], Y = new float[n] };
        for (int i = 0; i < n; i++)
        {
            g.X[i] = (float)((to[i] + from[i]) / 2);
            g.Y[i] = codes[i];
        }

        g.SetTruescale(g.X[0], g.X[n - 1], 0f, 1.37f);
        g.Xaxis = "Sequence mid position";
        g.Yaxis = "TESTCODE value";
        g.AddLine(g.X[0], 0.74f, g.X[n - 1], 0.74f, 1);
        g.AddLine(g.X[0], 0.95f, g.X[n - 1], 0.95f, 3);
        graph.Title = graph.InitialTitle;
        graph.AppendTitle(EmbossFeatureReports.MakeUsa(seq));
        graph.Replace(g);
        graph.Display(ctx);
    }

    private static EmbossReport TcodeReport(EmbossSeq seq, List<int> from, List<int> to, List<float> codes, bool firstTable)
    {
        var features = new List<EmbossFeature>();
        for (int i = 0; i < codes.Count; i++)
        {
            float fval = codes[i];
            string coding = (double)fval < 0.74 ? "*Estimation Non-coding" : (double)fval >= 0.95 ? "*Estimation Coding" : "*Estimation No opinion";
            features.Add(new EmbossFeature("misc_feature", "tcode", from[i], to[i], fval, '\0', 0, [new EmbossFeatureTag("note", coding)]));
        }

        var sequence = new EmbossReportSequence(seq.Name, EmbossFeatureReports.MakeUsa(seq), seq.Acc, seq.Desc, seq.Seq.ToString(), EmbossSeqRange.GetBegin(seq), EmbossSeqRange.GetEnd(seq), seq.Offset, seq.Reversed, seq.Type == "N");
        // the feature table is made for the first sequence (ajFeattableNewSeq) and cleared after every sequence (ajFeattableClear: no name, no type, no length: ajFeatNew raises the length to the largest end)
        if (firstTable)
            return new EmbossReport(sequence, features, false, "N", TableLength: sequence.Length, Tags: [new("str", "Estimation", "Estimation")]);
        long length = features.Count == 0 ? 0 : features.Max(f => (long)Math.Max(f.Start, f.End));
        return new EmbossReport(sequence, features, false, "", TableLength: length, Tags: [new("str", "Estimation", "Estimation")]) { TableSeqId = "", TypeUnknown = true };
    }
}
