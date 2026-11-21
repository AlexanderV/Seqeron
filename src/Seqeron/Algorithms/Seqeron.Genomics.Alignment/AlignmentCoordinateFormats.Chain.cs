using System.Text;
using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentCoordinateFormats
{
    /// <summary><c>Bio/Align/chain.py</c> <c>AlignmentWriter.format_alignment</c>.</summary>
    private static string FormatChain(BioCoordAlignment alignment)
    {
        var coordinates = Copy(alignment.Coordinates);
        if (coordinates.Length == 0 || coordinates[0].Length == 0)
            return "";
        var (target, query) = Pair(alignment);
        string qName = Id(query, "None"), tName = Id(target, "None");
        long tSize = target.Length, qSize = query.Length;
        var lines = new List<string>();
        string score = TryAttr(alignment, "score", out var sc) ? G(sc) : "0";
        string tStrand = "+", qStrand = "+";
        if (coordinates[0][0] > coordinates[0][^1])
        {
            tStrand = "-";
            for (int i = 0; i < coordinates[0].Length; i++)
                coordinates[0][i] = tSize - coordinates[0][i];
        }

        if (coordinates[1][0] > coordinates[1][^1])
        {
            qStrand = "-";
            for (int i = 0; i < coordinates[1].Length; i++)
                coordinates[1][i] = qSize - coordinates[1][i];
        }

        long tStart = coordinates[0][0], tEnd = coordinates[0][^1], qStart = coordinates[1][0], qEnd = coordinates[1][^1];
        string line = $"chain {score} {tName} {Str(tSize)} {tStrand} {Str(tStart)} {Str(tEnd)} {qName} {Str(qSize)} {qStrand} {Str(qStart)} {Str(qEnd)}";
        if (alignment.Attributes.TryGetPropertyValue("annotations", out var ann) && ann is JsonObject ao && ao.TryGetPropertyValue("id", out var cid) && cid is not null)
            line += " " + PyStr(cid);
        lines.Add(line);
        tStart = coordinates[0][0];
        qStart = coordinates[1][0];
        long step = 0, qGap = 0, tGap = 0;
        for (int k = 1; k < coordinates[0].Length; k++)
        {
            tEnd = coordinates[0][k];
            qEnd = coordinates[1][k];
            if (tStart == tEnd)
            {
                if (tGap > 0)
                {
                    qGap = qEnd - qStart;
                    lines.Add($"{Str(step)}\t{Str(tGap)}\t{Str(qGap)}");
                    step = 0;
                    tGap = 0;
                    qGap = 0;
                }
                else if (qGap > 0)
                {
                    lines.Add($"{Str(step)}\t{Str(tGap)}\t{Str(qGap)}");
                    step = 0;
                    qGap = qEnd - qStart;
                }
                else
                {
                    qGap = qEnd - qStart;
                }

                qStart = qEnd;
            }
            else if (qStart == qEnd)
            {
                if (qGap > 0)
                {
                    tGap = tEnd - tStart;
                    lines.Add($"{Str(step)}\t{Str(tGap)}\t{Str(qGap)}");
                    step = 0;
                    tGap = 0;
                    qGap = 0;
                }
                else if (tGap > 0)
                {
                    lines.Add($"{Str(step)}\t{Str(tGap)}\t{Str(qGap)}");
                    step = 0;
                    tGap = tEnd - tStart;
                }
                else
                {
                    tGap = tEnd - tStart;
                }

                tStart = tEnd;
            }
            else
            {
                long tStep, qStep;
                if (step == 0)
                {
                    tStep = tEnd - tStart;
                    qStep = qEnd - qStart;
                    if (tGap == 0 && qGap == 0)
                    {
                        step = tStep;
                    }
                    else
                    {
                        lines.Add($"{Str(tStep)}\t{Str(tGap)}\t{Str(qGap)}");
                        tGap = 0;
                        qGap = 0;
                    }
                }
                else
                {
                    lines.Add($"{Str(step)}\t{Str(tGap)}\t{Str(qGap)}");
                    tStep = tEnd - tStart;
                    qStep = qEnd - qStart;
                    step = tStep;
                    tGap = 0;
                    qGap = 0;
                }

                if (tStep != qStep)
                    throw new ArgumentException($"Expected equal step size in target and query (found {tStep} and {qStep})");
                tStart = tEnd;
                qStart = qEnd;
            }
        }

        if (tGap > 0 || qGap > 0)
            lines.Add($"0\t{Str(tGap)}\t{Str(qGap)}");
        lines.Add(Str(step));
        lines.Add("");
        return string.Join('\n', lines) + "\n";
    }

    /// <summary><c>Bio/Align/chain.py</c> <c>AlignmentIterator</c>.</summary>
    private static List<BioCoordAlignment> ParseChain(string text)
    {
        var result = new List<BioCoordAlignment>();
        var all = Lines(text).ToList();
        int p = 0;
        while (p < all.Count && all[p].StartsWith('#'))
            p++;
        string pending = p < all.Count ? all[p++] : "";
        while (pending.Length > 0)
        {
            var words = Words(pending);
            pending = "";
            string? chainId = null;
            if (words.Length == 13)
                chainId = words[12];
            else if (words.Length != 12)
                throw new ArgumentException("expected line with 12 or 13 columns at the beginning of a chain block");
            if (words[0] != "chain")
                throw new ArgumentException("expected first line of chain block to start with 'chain '");
            double score = ParseDouble(words[1]);
            string tName = words[2];
            long tSize = PyInt(words[3]);
            string tStrand = words[4];
            long tStart = PyInt(words[5]);
            string qName = words[7];
            long qSize = PyInt(words[8]);
            string qStrand = words[9];
            long qStart = PyInt(words[10]);
            long qPosition = 0, tPosition = 0;
            var qStarts = new List<long> { qPosition };
            var tStarts = new List<long> { tPosition };
            while (p < all.Count)
            {
                string line = all[p++];
                var w = Words(line);
                if (w.Length == 0)
                    throw new ArgumentException("list index out of range");
                if (w[0] == "chain")
                {
                    pending = line;
                    break;
                }

                long step = PyInt(w[0]);
                qPosition += step;
                tPosition += step;
                qStarts.Add(qPosition);
                tStarts.Add(tPosition);
                if (w.Length == 1)
                    break;
                long tGap = PyInt(w[1]);
                long qGap = PyInt(w[2]);
                if (tGap > 0)
                {
                    tPosition += tGap;
                    qStarts.Add(qPosition);
                    tStarts.Add(tPosition);
                }

                if (qGap > 0)
                {
                    qPosition += qGap;
                    qStarts.Add(qPosition);
                    tStarts.Add(tPosition);
                }
            }

            var co = new[] { tStarts.ToArray(), qStarts.ToArray() };
            for (int i = 0; i < co[0].Length; i++)
            {
                co[0][i] += tStart;
                co[1][i] += qStart;
            }

            if (tStrand == "-")
                for (int i = 0; i < co[0].Length; i++) co[0][i] = tSize - co[0][i];
            else if (tStrand != "+")
                throw new ArgumentException($"tStrand should be '+' or '-' (found '{tStrand}')");
            if (qStrand == "-")
                for (int i = 0; i < co[1].Length; i++) co[1][i] = qSize - co[1][i];
            else if (qStrand != "+")
                throw new ArgumentException($"qStrand should be '+' or '-' (found '{qStrand}')");
            var attrs = new JsonObject { ["score"] = score };
            if (chainId != null)
                attrs["annotations"] = new JsonObject { ["id"] = chainId };
            result.Add(new BioCoordAlignment { Sequences = [Unknown(tName, tSize), Unknown(qName, qSize)], Coordinates = co, Attributes = attrs });
            if (pending.Length == 0)
            {
                // a blank line between the blocks is optional
                string line = p < all.Count ? all[p++] : "";
                if (line.Trim().Length == 0)
                    line = p < all.Count ? all[p++] : "";
                pending = line;
            }
        }

        return result;
    }

    private static double ParseDouble(string s)
    {
        if (!double.TryParse(s.Trim(), System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out double v))
            throw new ArgumentException($"could not convert string to float: '{s}'");
        return v;
    }
}
