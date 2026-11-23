using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentOperations
{
    // the residues [start, end) of a record, '?' where they are undefined (Python: bytes(seq) after UndefinedSequenceError, filled from defined_ranges)
    private static string DefinedText(BioCoordSequence s, long start, long end)
    {
        var chars = new char[Math.Max(0, end - start)];
        Array.Fill(chars, '?');
        if (s.Parts is not null)
        {
            foreach (var (partStart, text) in s.Parts)
            {
                long a = Math.Max(start, partStart), b = Math.Min(end, partStart + text.Length);
                for (long p = a; p < b; p++)
                    chars[p - start] = text[(int)(p - partStart)];
            }
        }

        return new string(chars);
    }

    private static string PositionText(string name, long start, int positionWidth)
    {
        string positionText = start.ToString(CultureInfo.InvariantCulture);
        int offset = positionWidth - positionText.Length - 1;
        if (offset < 0)
            return name + " .." + positionText[(-offset + 3)..] + " ";
        return name + new string(' ', offset) + positionText + " ";
    }

    private static string Pad9(long v) => v.ToString(CultureInfo.InvariantCulture).PadLeft(9);

    /// <summary>
    /// <c>str(alignment)</c> / <c>alignment.format()</c>: Biopython's "pretty" text of an alignment: for two sequences the target and query rows with the match line (<c>|</c> identical, <c>.</c> different, <c>-</c> gap,
    /// <c>:</c> positive in <paramref name="matrix"/>), 80 columns wide with the sequence positions at the line ends; for more sequences the rows only; a protein aligned to nucleotides (coordinates in steps of 1 and 3) is written
    /// with the codons spread over three columns; a record's identifier (first 9 characters) names the row, <c>target</c> / <c>query</c> otherwise; undefined residues are <c>?</c>.
    /// </summary>
    /// <param name="alignment">The alignment.</param>
    /// <param name="matrix">Optional substitution matrix for the <c>:</c> positions.</param>
    public static string Format(BioCoordAlignment alignment, SubstitutionMatrix? matrix = null)
    {
        RequireCoordinates(alignment);
        int n = alignment.Sequences.Count;
        bool writePattern = n == 2;
        var coordinates = alignment.Coordinates;
        var steps = StepsOf(coordinates);
        var aligned = AlignedBlocks(steps);
        const int nameWidth = 10;
        var names = new List<string>();
        var seqs = new List<string>();
        var indices = coordinates.Select(r => new long[r.Length]).ToArray();
        for (int i = 0; i < n; i++)
        {
            var record = alignment.Sequences[i];
            string name;
            if (record.Id is not null)
                name = record.Id.Length > nameWidth - 1 ? record.Id[..(nameWidth - 1)] : record.Id;
            else if (n == 2)
                name = i == 0 ? "target" : "query";
            else
                name = "";
            names.Add(name.PadRight(nameWidth));
            long[] positions = coordinates[i];
            long start = positions.Min(), end = positions.Max();
            string seq = DefinedText(record, start, end);
            long[] alignedSteps = [.. steps[i].Where((_, k) => aligned[k])];
            if (alignedSteps.Length == 0)
                alignedSteps = steps[i];
            if (alignedSteps.Count(x => x > 0) >= alignedSteps.Count(x => x < 0))
            {
                for (int k = 0; k < positions.Length; k++)
                    indices[i][k] = positions[k] - start;
            }
            else
            {
                for (int k = 0; k < steps[i].Length; k++)
                    steps[i][k] = -steps[i][k];
                seq = AlignmentCoordinateFormats.ReverseComplement(seq);
                end = positions.Max();
                for (int k = 0; k < positions.Length; k++)
                    indices[i][k] = end - positions[k];
            }

            if (seq.Any(c => c > 127))
                return FormatUnicode(alignment, matrix);
            seqs.Add(seq);
        }

        int blocks = steps.Length == 0 ? 0 : steps[0].Length;
        var merged = new long[blocks];
        for (int k = 0; k < blocks; k++)
        {
            long min = long.MaxValue, max = long.MinValue;
            for (int i = 0; i < n; i++)
            {
                min = Math.Min(min, steps[i][k]);
                max = Math.Max(max, steps[i][k]);
            }

            merged[k] = -min > max ? min : max;
        }

        for (int i = 0; i < n; i++)
        {
            var row = indices[i];
            var rowSteps = new List<long>();
            var alignedSteps = new List<long>();
            for (int k = 0; k < blocks; k++)
            {
                long rs = row[k + 1] - row[k];
                if (rs > 0 && aligned[k])
                {
                    rowSteps.Add(rs);
                    alignedSteps.Add(merged[k]);
                }
            }

            if (rowSteps.Zip(alignedSteps).All(p => p.First == p.Second))
            {
                // one letter per column
            }
            else if (rowSteps.Zip(alignedSteps).All(p => 3 * p.First == p.Second))
            {
                for (int k = 0; k < row.Length; k++)
                    row[k] *= 3;
                seqs[i] = string.Join("  ", seqs[i].Select(c => c.ToString())) + "  ";
                writePattern = false;
            }
            else
            {
                throw Bad("Inconsistent coordinates");
            }
        }

        const int prefixWidth = 10, positionWidth = 10, lineWidth = 80;
        var lines = new List<string>();
        var idxSteps = new long[blocks];
        for (int k = 0; k < blocks; k++)
        {
            long min = long.MaxValue, max = long.MinValue;
            for (int i = 0; i < n; i++)
            {
                long d = indices[i][k + 1] - indices[i][k];
                min = Math.Min(min, d);
                max = Math.Max(max, d);
            }

            idxSteps[k] = -min > max ? min : max;
        }

        int column = lineWidth;
        for (int i = 0; i < n; i++)
        {
            string name = names[i], seq = seqs[i];
            long[] positions = coordinates[i], row = indices[i];
            long start = positions[0];
            column = lineWidth;
            long startIndex = row[0];
            for (int k = 0; k < blocks; k++)
            {
                long step = idxSteps[k], end = positions[k + 1], endIndex = row[k + 1];
                if (step < 0)
                {
                    if (prefixWidth + positionWidth < column)
                    {
                        string positionText = start.ToString(CultureInfo.InvariantCulture);
                        int offset = positionWidth - positionText.Length - 1;
                        lines[^1] += offset < 0 ? " .." + positionText[(-offset + 3)..] : " " + positionText;
                    }

                    column = lineWidth;
                    start = end;
                    startIndex = endIndex;
                    continue;
                }

                string s = endIndex == startIndex ? new string('-', (int)step) : Sub(seq, startIndex, endIndex);
                while (column + s.Length >= lineWidth)
                {
                    int rest = lineWidth - column;
                    if (rest > 0)
                    {
                        lines[^1] += s[..Math.Min(rest, s.Length)];
                        s = s[Math.Min(rest, s.Length)..];
                        if (start != end)
                        {
                            long advance = endIndex - startIndex == Math.Abs(end - start) ? rest : (rest + 2) / 3;
                            start = start < end ? start + advance : start - advance;
                        }

                        startIndex += rest;
                    }

                    lines.Add(PositionText(name, start, positionWidth));
                    column = nameWidth + positionWidth;
                }

                lines[^1] += s;
                if (startIndex != endIndex)
                {
                    startIndex = endIndex;
                    start = end;
                }

                column += s.Length;
            }
        }

        long[] finals = [.. coordinates.Select(r => r[^1])];
        if (writePattern)
        {
            int m = lines.Count / 2;
            var lines1 = lines.Take(m).ToList();
            var lines2 = lines.Skip(m).ToList();
            var patternLines = new List<string>();
            long position = 0;
            for (int p = 0; p < Math.Min(lines1.Count, lines2.Count); p++)
            {
                string aligned1 = lines1[p][(nameWidth + positionWidth)..], aligned2 = lines2[p][(nameWidth + positionWidth)..];
                var pattern = new StringBuilder();
                for (int c = 0; c < Math.Min(aligned1.Length, aligned2.Length); c++)
                {
                    char c1 = aligned1[c], c2 = aligned2[c];
                    char symbol;
                    if (c1 == c2)
                    {
                        if (c1 == ' ')
                            break;
                        symbol = '|';
                    }
                    else if (c1 == '-' || c2 == '-')
                    {
                        symbol = '-';
                    }
                    else
                    {
                        symbol = '.';
                        if (matrix is not null && c1 != ' ' && c2 != ' ')
                        {
                            if (!matrix.TryGetValue(char.ToUpperInvariant(c1), char.ToUpperInvariant(c2), out decimal value))
                                throw Bad("index out of range: letter not in the substitution matrix");
                            if (value > 0)
                                symbol = ':';
                        }
                    }

                    pattern.Append(symbol);
                }

                patternLines.Add("          " + Pad9(position) + " " + pattern);
                position += pattern.Length;
            }

            int finalPositionWidth = Math.Max(finals.Max(), position).ToString(CultureInfo.InvariantCulture).Length;
            if (column + finalPositionWidth <= lineWidth)
            {
                if (prefixWidth + positionWidth < column)
                {
                    lines1[^1] += " " + finals[0].ToString(CultureInfo.InvariantCulture).PadLeft(finalPositionWidth);
                    lines2[^1] += " " + finals[1].ToString(CultureInfo.InvariantCulture).PadLeft(finalPositionWidth);
                    patternLines[^1] += " " + position.ToString(CultureInfo.InvariantCulture).PadLeft(finalPositionWidth);
                }
            }
            else
            {
                lines1.Add(names[0] + finals[0].ToString(CultureInfo.InvariantCulture).PadLeft(9));
                patternLines.Add("          " + Pad9(position));
                lines2.Add(names[1] + finals[1].ToString(CultureInfo.InvariantCulture).PadLeft(9));
            }

            var result = new List<string>();
            for (int p = 0; p < Math.Min(lines1.Count, Math.Min(lines2.Count, patternLines.Count)); p++)
                result.Add($"{lines1[p]}\n{patternLines[p]}\n{lines2[p]}\n");
            return string.Join("\n", result);
        }
        else
        {
            int m = lines.Count / n;
            int finalPositionWidth = finals.Max().ToString(CultureInfo.InvariantCulture).Length;
            var blockTexts = new List<string>();
            if (column + finalPositionWidth < lineWidth)
            {
                if (prefixWidth + positionWidth < column)
                {
                    for (int i = 0; i < n; i++)
                        lines[m - 1 + i * m] += " " + finals[i].ToString(CultureInfo.InvariantCulture).PadLeft(finalPositionWidth);
                }

                for (int j = 0; j < m; j++)
                    blockTexts.Add(string.Join("\n", Enumerable.Range(0, n).Select(i => lines[j + i * m])) + "\n");
            }
            else
            {
                for (int j = 0; j < m; j++)
                    blockTexts.Add(string.Join("\n", Enumerable.Range(0, n).Select(i => lines[j + i * m])) + "\n");
                blockTexts.Add(string.Join("\n", Enumerable.Range(0, n).Select(i => names[i] + finals[i].ToString(CultureInfo.InvariantCulture).PadLeft(9))) + "\n");
            }

            return string.Join("\n", blockTexts);
        }
    }

    // Alignment._format_unicode: the plain aligned rows (non-ASCII residues)
    private static string FormatUnicode(BioCoordAlignment alignment, SubstitutionMatrix? matrix)
    {
        var seqs = new List<string>();
        var coordinates = alignment.Coordinates.Select(r => (long[])r.Clone()).ToArray();
        for (int i = 0; i < coordinates.Length; i++)
        {
            string seq = Contents(alignment.Sequences[i]);
            long[] row = coordinates[i];
            if (row[0] > row[^1])
            {
                for (int k = 0; k < row.Length; k++)
                    row[k] = seq.Length - row[k];
                seq = AlignmentCoordinateFormats.ReverseComplement(seq);
            }

            seqs.Add(seq);
        }

        var steps = StepsOf(coordinates);
        int blocks = steps.Length == 0 ? 0 : steps[0].Length;
        var maxSteps = new long[blocks];
        for (int k = 0; k < blocks; k++)
            maxSteps[k] = steps.Max(r => r[k]);
        var alignedSeqs = new List<string>();
        for (int i = 0; i < coordinates.Length; i++)
        {
            var sb = new StringBuilder();
            long start = coordinates[i][0];
            for (int k = 0; k < blocks; k++)
            {
                long end = coordinates[i][k + 1];
                sb.Append(end == start ? new string('-', (int)maxSteps[k]) : Sub(seqs[i], start, end));
                start = end;
            }

            alignedSeqs.Add(sb.ToString());
        }

        if (seqs.Count > 2)
            return string.Join("\n", alignedSeqs) + "\n";
        string a = alignedSeqs[0], b = alignedSeqs[1];
        var pattern = new StringBuilder();
        for (int c = 0; c < Math.Min(a.Length, b.Length); c++)
        {
            char c1 = a[c], c2 = b[c];
            char symbol;
            if (c1 == c2)
            {
                symbol = '|';
            }
            else if (c1 == '-' || c2 == '-')
            {
                symbol = '-';
            }
            else
            {
                symbol = '.';
                if (matrix is not null && c1 != ' ' && c2 != ' ' && matrix.TryGetValue(char.ToUpperInvariant(c1), char.ToUpperInvariant(c2), out decimal value) && value > 0)
                    symbol = ':';
            }

            pattern.Append(symbol);
        }

        return $"{a}\n{pattern}\n{b}\n";
    }
}
