using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// <c>nucmer2xfig</c> of MUMmer 3.23 (the perl script of Steven Salzberg, 2001): an xfig plot of a coords file made by <c>show-coords -c -l</c> (the contigs of one genome against a reference chromosome): the reference as
/// a red horizontal line with coordinate and percent identity tick marks, every aligned stretch as a green line at the height of its identity, every contig (the stretch its alignments span, split where the alignments are
/// more than twice the contig length apart) as a blue line with a vertical label when the line is wider than the label.
/// </summary>
public static class MummerNucmer2Xfig
{
    private const double XScale = 0.005, YScale = 20;
    private const int ChrColor = 4, Green = 2, ContigColor = 1;

    /// <summary>The numeric value perl gives a string (<c>$x + 0</c>): the leading number, 0 when there is none.</summary>
    private static double Num(string? s)
    {
        if (string.IsNullOrEmpty(s))
            return 0;
        int i = 0;
        while (i < s.Length && char.IsWhiteSpace(s[i]))
            i++;
        int start = i;
        if (i < s.Length && (s[i] == '+' || s[i] == '-'))
            i++;
        int digits = i;
        while (i < s.Length && char.IsAsciiDigit(s[i]))
            i++;
        if (i < s.Length && s[i] == '.')
        {
            i++;
            while (i < s.Length && char.IsAsciiDigit(s[i]))
                i++;
        }

        if (i == digits || (i == digits + 1 && s[digits] == '.'))
            return 0;
        if (i < s.Length && (s[i] == 'e' || s[i] == 'E'))
        {
            int j = i + 1;
            if (j < s.Length && (s[j] == '+' || s[j] == '-'))
                j++;
            int e = j;
            while (j < s.Length && char.IsAsciiDigit(s[j]))
                j++;
            if (j > e)
                i = j;
        }

        return double.Parse(s.AsSpan(start, i - start), NumberStyles.Float, CultureInfo.InvariantCulture);
    }

    /// <summary>perl's <c>int ()</c>: truncation toward zero.</summary>
    private static double Int(double v) => Math.Truncate(v);

    private static string F0(double v) => CFormat.Fixed(v, 0);

    /// <summary>perl's string form of a number (<c>"$x"</c>): integers without a decimal point, otherwise 15 significant digits.</summary>
    private static string Str(double v) => v == Math.Truncate(v) && Math.Abs(v) < 1e15 ? ((long)v).ToString(CultureInfo.InvariantCulture) : CFormat.General(v, 15);

    /// <summary>Runs <c>nucmer2xfig coordsfile</c>.</summary>
    /// <param name="coords">The text of the coords file (<c>show-coords -c -l</c> without its header, one line per alignment of 19 white-space separated fields).</param>
    public static string Run(string coords)
    {
        ArgumentNullException.ThrowIfNull(coords);
        var sb = new StringBuilder("#FIG 3.2\nLandscape\nCenter\nInches\nLetter  \n100.00\nSingle\n-2\n1200 2\n");
        var left = new Dictionary<string, double>();
        var right = new Dictionary<string, double>();
        bool firstTime = true;
        string prevQid = string.Empty;
        double prevS1 = 0, s1 = 0, e1 = 0, rlen = 0, qlen = 0, percent = 0;
        string qid = string.Empty, rid = string.Empty;

        var lines = coords.Split('\n').ToList();
        if (lines.Count > 0 && lines[^1].Length == 0)
            lines.RemoveAt(lines.Count - 1);
        foreach (string raw in lines)
        {
            var f = raw.Split([' ', '\t', '\n', '\r', '\f', '\v'], StringSplitOptions.RemoveEmptyEntries);
            string G(int i) => i < f.Length ? f[i] : string.Empty;
            s1 = Num(G(0));
            e1 = Num(G(1));
            percent = Num(G(9));
            rlen = Num(G(11));
            qlen = Num(G(12));
            rid = G(17);
            qid = G(18);
            if (prevQid == qid)
            {
                double dist = Math.Abs(s1 - prevS1);
                if (dist > 2 * qlen)
                {
                    double l = left.GetValueOrDefault(qid), r = right.GetValueOrDefault(qid);
                    if (r - l > 1000)
                    {
                        Line(sb, l, 50, r, 50, ContigColor);
                        Label(sb, l, 50, r, qid, 5);
                    }

                    left[qid] = s1;
                    right[qid] = e1;
                }
                else
                {
                    if (s1 < left.GetValueOrDefault(qid))
                        left[qid] = s1;
                    if (e1 > right.GetValueOrDefault(qid))
                        right[qid] = e1;
                }
            }
            else
            {
                left[qid] = s1;
                right[qid] = e1;
                if (!firstTime)
                {
                    Line(sb, left.GetValueOrDefault(prevQid), 50, right.GetValueOrDefault(prevQid), 50, ContigColor);
                    Label(sb, left.GetValueOrDefault(prevQid), 50, right.GetValueOrDefault(prevQid), prevQid, 5);
                }
                else
                    firstTime = false;
            }

            prevQid = qid;
            prevS1 = s1;
            if (percent < 50)
                percent = 50;
            Line(sb, s1, percent - 50, e1, percent - 50, Green);
        }

        left[qid] = s1;
        right[qid] = e1;
        Line(sb, s1, 50, e1, 50, ContigColor);
        Label(sb, s1, 50, e1, qid, 5);

        double labelX = XScale * rlen;
        sb.Append("4 0 0 100 0 0 12 0.0000 4 135 405 ").Append(F0(labelX)).Append(' ').Append(F0(0)).Append(' ').Append(rid).Append("\\001\n");
        sb.Append("2 1 0 2 ").Append(ChrColor).Append(" 7 100 0 -1 0.000 0 0 -1 0 0 2\n");
        sb.Append("\t ").Append(F0(0)).Append(' ').Append(F0(0)).Append(' ').Append(F0(XScale * rlen)).Append(' ').Append(F0(0)).Append('\n');
        const int pointSize = 5;
        for (double i = 250000; i < rlen - 50000; i += 250000)
        {
            sb.Append("4 0 0 100 0 0 ").Append(pointSize).Append(" 0.0000 4 135 405 ");
            sb.Append(F0((i * XScale) + 20)).Append(' ').Append(F0(-20 + (YScale * -0.5)));
            sb.Append(' ').Append(Str(i)).Append("\\001\n");
            sb.Append("2 1 0 1 0 7 50 0 -1 0.000 0 0 -1 0 0 2\n");
            sb.Append(F0(i * XScale)).Append(' ').Append(F0(0)).Append(' ').Append(F0(i * XScale)).Append(" -50\n");
        }

        for (int p = 50; p < 101; p += 10)
        {
            sb.Append("4 0 0 100 0 0 ").Append(pointSize).Append(" 0.0000 4 135 405 ");
            sb.Append("-150 ").Append(F0(((p - 50) * YScale) + 10));
            sb.Append(' ').Append(p.ToString(CultureInfo.InvariantCulture)).Append("\\001\n");
            sb.Append("2 1 0 1 0 7 50 0 -1 0.000 0 0 -1 0 0 2\n");
            sb.Append("-50 ").Append(F0((p - 50) * YScale)).Append(" 0 ").Append(F0((p - 50) * YScale)).Append('\n');
        }

        return sb.ToString();
    }

    private static void Line(StringBuilder sb, double xleft, double yleft, double xright, double yright, int color)
    {
        double xl = Int(XScale * xleft), xr = Int(XScale * xright);
        if (xl == xr)
            xr += 1;
        sb.Append("2 1 0 2 ").Append(color).Append(" 7 100 0 -1 0.000 0 0 -1 0 0 2\n");
        sb.Append("\t ").Append(F0(xl)).Append(' ').Append(F0(YScale * yleft)).Append(' ').Append(F0(xr)).Append(' ').Append(F0(YScale * yright)).Append('\n');
    }

    private static void Label(StringBuilder sb, double xleft, double yleft, double xright, string id, int psize)
    {
        double y = (yleft * YScale) + 20;
        double xl = xleft * XScale, xr = xright * XScale;
        double textHeight = psize * 8;
        if (xr - xl > textHeight)
        {
            sb.Append("4 0 0 50 0 0 ").Append(psize).Append(" 4.7124 4 135 435 ");
            sb.Append(F0(xl)).Append(' ').Append(F0(y)).Append(' ').Append(id).Append("\\001\n");
        }
    }
}
