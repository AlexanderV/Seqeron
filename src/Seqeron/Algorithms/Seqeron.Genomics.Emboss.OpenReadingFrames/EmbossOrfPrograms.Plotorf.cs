// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.OpenReadingFrames, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127, S1905
using System.Globalization;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.OpenReadingFrames;

using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Emboss.SequenceFormats.EmbossSequencePrograms;

/// <summary><c>plotorf</c> of EMBOSS 6.6.0 (<c>emboss/plotorf.c</c>): the potential open reading frames of a nucleotide sequence as a plot of six frames (B09 WP63, F282).</summary>
public static partial class EmbossOrfPrograms
{
    // ajArrCommaList: n = ajStrParseCountC(s, ",\n") counts the tokens from the first character that is no delimiter (ajStrTokenAssignC skips the leading delimiters), but the tokens are parsed from the start (ajStrTokenNewC):
    // a leading delimiter gives an empty first token and the last token(s) beyond n are lost; each token is trimmed of white space
    internal static List<string> CommaList(string s)
    {
        int n = 0;
        int pos = 0;
        while (pos < s.Length && s[pos] is ',' or '\n')
            pos++;
        while (pos < s.Length)
        {
            n++;
            while (pos < s.Length && s[pos] is not (',' or '\n'))
                pos++;
            while (pos < s.Length && s[pos] is ',' or '\n')
                pos++;
        }

        var tokens = new List<string>();
        pos = 0;
        while (pos < s.Length && tokens.Count < n)
        {
            int end = s.IndexOfAny([',', '\n'], pos);
            if (end < 0)
                end = s.Length;
            tokens.Add(s[pos..end].Trim(' ', '\t', '\n', '\v', '\f', '\r'));
            pos = end;
            while (pos < s.Length && s[pos] is ',' or '\n')
                pos++;
        }

        return tokens;
    }

    // strncmp(p, codon, 3) == 0 for every codon of the list: up to three characters, equal ones and the NUL that ends both strings
    private static bool PlotorfIsIn(string p, int at, List<string> codons)
    {
        foreach (string c in codons)
        {
            bool equal = true;
            for (int k = 0; k < 3 && equal; k++)
            {
                char x = at + k < p.Length ? p[at + k] : '\0';
                char y = k < c.Length ? c[k] : '\0';
                equal = x == y;
                if (x == '\0')
                    break;
            }

            if (equal)
                return true;
        }

        return false;
    }

    // plotorf_norfs: the start and end of every potential ORF of the frame n
    private static (float[] X, float[] Y) PlotorfNorfs(string seq, string rev, int n, int beg, List<string> starts, List<string> stops)
    {
        int len = seq.Length;
        string p;
        int po;
        if (n < 3)
        {
            p = seq;
            po = n % 3;
        }
        else
        {
            p = rev;
            po = len % 3;
            po -= n % 3;
            if (po < 0)
                po += 3;
        }

        bool inframe = false;
        int count = 0;
        for (int i = po; i < len - 2; i += 3)
        {
            if (PlotorfIsIn(p, i, starts) && !inframe)
            {
                ++count;
                inframe = true;
                continue;
            }

            if (PlotorfIsIn(p, i, stops) && inframe)
                inframe = false;
        }

        var x = new float[count];
        var y = new float[count];
        count = 0;
        inframe = false;
        for (int i = po; i < len - 2; i += 3)
        {
            if (PlotorfIsIn(p, i, starts) && !inframe)
            {
                if (x.Length > 0)
                    x[count] = n < 3 ? i + beg : len - i - 1 + beg;
                ++count;
                inframe = true;
                continue;
            }

            if (PlotorfIsIn(p, i, stops) && inframe)
            {
                if (x.Length > 0)
                    y[count - 1] = n < 3 ? i + beg : len - i - 1 + beg;
                inframe = false;
            }
        }

        if (inframe && x.Length > 0)
            y[count - 1] = n < 3 ? len + beg - 1 : beg;
        return (x, y);
    }

    /// <summary><c>plotorf</c>: for each of the six frames the stretches from a start codon to the next stop codon of the frame (the codons of <c>-start</c> and <c>-stop</c>), as the rectangles of six plots (<c>-graph data</c> writes them as six data files).</summary>
    internal static EmbossRunResult Plotorf(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("plotorf");
        var all = EmbossProgramBase.Open("plotorf", a, files, ctx, "dna", single: true);
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        var seq = all.Peek();
        EmbossGraphXy? graph;
        string start, stop;
        try
        {
            start = a.String("start", "ATG");
            stop = a.String("stop", "TAA,TAG,TGA");
            graph = EmbossGraphArgs.Open(a, ctx, "plotorf", "graph", "", 6, EmbossOutPath.DefaultBaseName(seq.Name, a.Raw("sequence") ?? ""), files);
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        if (graph is null)
            return EmbossProgramBase.Died(ctx, "graph");
        var starts = CommaList(start.ToUpperInvariant());
        var stops = CommaList(stop.ToUpperInvariant());
        int beg = EmbossSeqRange.GetBegin(seq), end = EmbossSeqRange.GetEnd(seq);
        string upper = seq.SeqText.ToUpperInvariant();
        string str = EmbossProgramBase.SubS(upper, beg - 1, end - 1);
        string rev = EmbossNuc.ReverseComplement(str);
        string[] titles = ["F1", "F2", "F3", "R1", "R2", "R3"];
        for (int i = 0; i < 6; ++i)
        {
            var (x, y) = PlotorfNorfs(str, rev, i, beg, starts, stops);
            var data = new EmbossGraphData { Type = "Multi 2D Plot Small", Title = titles[i], Xaxis = "Sequence", Yaxis = "Orf" };
            graph.Replace(data, i);
            data.SetTruescale(beg, end, 0f, 1f);
            for (int j = 0; j < x.Length; ++j)
                data.Objects.Add(new EmbossGraphObject("Filled Rectangle", y[j], 0f, x[j], 1f, 4));
        }

        foreach (var g in graph.Graphs)
        {
            g.MinX = beg;
            g.MaxX = end;
            g.MinY = 0f;
            g.MaxY = 1f;
        }

        graph.Title = "Potential codons (rectangles)";
        graph.Display(ctx, true);
        return new EmbossRunResult(graph.Stdout, EmbossRunResult.Join(ctx.Messages), 0, new Dictionary<string, string>(graph.Files));
    }
}
