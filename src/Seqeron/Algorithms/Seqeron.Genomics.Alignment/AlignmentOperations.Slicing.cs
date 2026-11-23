using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentOperations
{
    // numpy searchsorted on the (non-decreasing) cumulative gaps
    private static int SearchRight(long[] a, long value, int from = 0)
    {
        int i = from;
        while (i < a.Length && a[i] <= value)
            i++;
        return i;
    }

    private static int SearchLeft(long[] a, long value, int from = 0)
    {
        int i = from;
        while (i < a.Length && a[i] < value)
            i++;
        return i;
    }

    private static long[] CumulativeSum(long[] gaps)
    {
        var result = new long[gaps.Length];
        long sum = 0;
        for (int i = 0; i < gaps.Length; i++)
            result[i] = sum += gaps[i];
        return result;
    }

    private static int NormalizeColumn(int col, int m)
    {
        if (col < 0)
            col += m;
        if (col < 0 || col >= m)
            throw new ArgumentException($"column index {col} is out of bounds ({m} columns)");
        return col;
    }

    /// <summary><c>alignment[row, col]</c> with two integers: the letter (or <c>-</c>) of a row in a column.</summary>
    public static char Letter(BioCoordAlignment alignment, int row, int column)
    {
        var rows = Rows(alignment);
        row = NormalizeRow(row, rows.Length);
        column = NormalizeColumn(column, rows.Length == 0 ? 0 : rows[0].Length);
        return rows[row][column];
    }

    private static int NormalizeRow(int row, int n)
    {
        if (row < 0)
            row += n;
        if (row < 0 || row >= n)
            throw new ArgumentException("row index out of range");
        return row;
    }

    /// <summary><c>alignment[row, a:b:c]</c>: the printed row of one sequence restricted to the columns of a slice.</summary>
    public static string Row(BioCoordAlignment alignment, int row, PySlice columns)
    {
        var rows = Rows(alignment);
        row = NormalizeRow(row, rows.Length);
        string line = rows[row];
        var indices = columns.Range(line.Length);
        var sb = new StringBuilder();
        foreach (int i in indices)
            sb.Append(line[i]);
        return sb.ToString();
    }

    /// <summary><c>alignment[row, [c1, c2, ...]]</c>: the printed row of one sequence restricted to the listed columns (negative = from the end).</summary>
    public static string Row(BioCoordAlignment alignment, int row, IReadOnlyList<int> columns)
    {
        ArgumentNullException.ThrowIfNull(columns);
        var rows = Rows(alignment);
        row = NormalizeRow(row, rows.Length);
        string line = rows[row];
        var sb = new StringBuilder();
        foreach (int c in columns)
        {
            int i = c < 0 ? c + line.Length : c;
            if (i < 0 || i >= line.Length)
                throw new ArgumentException("string index out of range");
            sb.Append(line[i]);
        }

        return sb.ToString();
    }

    /// <summary><c>alignment[a:b:c, col]</c>: the letters (or <c>-</c>) of the selected rows in one column.</summary>
    public static string Column(BioCoordAlignment alignment, PySlice rows, int column)
    {
        var all = Rows(alignment);
        int m = all.Length == 0 ? 0 : all[0].Length;
        column = NormalizeColumn(column, m);
        var sb = new StringBuilder();
        foreach (int r in rows.Range(all.Length))
            sb.Append(all[r][column]);
        return sb.ToString();
    }

    /// <summary><c>alignment[a:b:c]</c>: the alignment of the selected sequences (the score is kept when the coordinates are unchanged).</summary>
    public static BioCoordAlignment Slice(BioCoordAlignment alignment, PySlice rows)
    {
        RequireCoordinates(alignment);
        var selected = rows.Range(alignment.Sequences.Count);
        var coordinates = selected.Select(i => (long[])alignment.Coordinates[i].Clone()).ToArray();
        var result = new BioCoordAlignment { Sequences = [.. selected.Select(i => alignment.Sequences[i])], Coordinates = coordinates };
        if (SameCoordinates(alignment.Coordinates, coordinates))
            return new BioCoordAlignment { Sequences = result.Sequences, Coordinates = coordinates, Attributes = ScoreOnly(alignment) };
        return result;
    }

    private static System.Text.Json.Nodes.JsonObject ScoreOnly(BioCoordAlignment a)
    {
        var attrs = new System.Text.Json.Nodes.JsonObject();
        if (a.Attributes.TryGetPropertyValue("score", out var s) && s is not null)
            attrs["score"] = s.DeepClone();
        return attrs;
    }

    private static bool SameCoordinates(long[][] a, long[][] b) => a.Length == b.Length && a.Zip(b).All(p => p.First.SequenceEqual(p.Second));

    /// <summary><c>alignment[a:b:c, d:e]</c>: the selected sequences restricted to the columns of a slice (a step of 1 and a non-empty range cut the blocks; other slices are taken as column lists).</summary>
    public static BioCoordAlignment Slice(BioCoordAlignment alignment, PySlice rows, PySlice columns)
    {
        var o = Orient(alignment);
        var gaps = Gaps(o.Steps);
        int m = (int)gaps.Sum();
        var (start, stop, step) = columns.Indices(m);
        if (start < stop && step == 1)
            return SliceRowsColumnRange(alignment, o, gaps, rows, start, stop);
        return SliceRowsColumnList(alignment, o, gaps, rows, columns.Range(m));
    }

    /// <summary><c>alignment[a:b:c, [c1, c2, ...]]</c>: the selected sequences restricted to a list of columns (the letters of the columns, re-aligned: the new sequences are the selected letters without gaps).</summary>
    public static BioCoordAlignment Slice(BioCoordAlignment alignment, PySlice rows, IReadOnlyList<int> columns)
    {
        ArgumentNullException.ThrowIfNull(columns);
        var o = Orient(alignment);
        var gaps = Gaps(o.Steps);
        return SliceRowsColumnList(alignment, o, gaps, rows, columns);
    }

    // _get_rows_cols_slice
    private static BioCoordAlignment SliceRowsColumnRange(BioCoordAlignment alignment, Oriented o, long[] gaps, PySlice rows, int startIndex, int stopIndex)
    {
        var selected = rows.Range(alignment.Sequences.Count);
        var coordinates = selected.Select(r => (long[])o.Coordinates[r].Clone()).ToArray();
        var steps = selected.Select(r => o.Steps[r]).ToArray();
        var rcs = selected.Select((r, k) => !coordinates[k].SequenceEqual(alignment.Coordinates[r])).ToArray();
        var indices = CumulativeSum(gaps);
        int i = SearchRight(indices, startIndex);
        int j = SearchLeft(indices, stopIndex, i) + 1;
        for (int r = 0; r < coordinates.Length; r++)
        {
            long offset = steps[r][i] - indices[i] + startIndex;
            if (steps[r][i] > 0)
                coordinates[r][i] += offset;
        }

        for (int r = 0; r < coordinates.Length; r++)
        {
            long offset = indices[j - 1] - stopIndex;
            if (steps[r][j - 1] > 0)
                coordinates[r][j] -= offset;
        }

        var cut = coordinates.Select(row => row[i..(j + 1)]).ToArray();
        var sequences = selected.Select(r => alignment.Sequences[r]).ToList();
        for (int r = 0; r < cut.Length; r++)
        {
            if (rcs[r])
                for (int k = 0; k < cut[r].Length; k++)
                    cut[r][k] = sequences[r].Length - cut[r][k];
        }

        var result = new BioCoordAlignment { Sequences = sequences, Coordinates = cut };
        if (SameCoordinates(alignment.Coordinates, cut))
            return new BioCoordAlignment { Sequences = sequences, Coordinates = cut, Attributes = ScoreOnly(alignment) };
        return result;
    }

    // _get_rows_cols_iterable
    private static BioCoordAlignment SliceRowsColumnList(BioCoordAlignment alignment, Oriented o, long[] gaps, PySlice rows, IReadOnlyList<int> columns)
    {
        var selected = rows.Range(alignment.Sequences.Count);
        var lines = new List<string>();
        foreach (int r in selected)
        {
            string text = Contents(o.Sequences[r]);
            var sb = new StringBuilder();
            long k = o.Coordinates[r][0];
            for (int c = 0; c < gaps.Length; c++)
            {
                long step = o.Steps[r][c];
                if (step != 0)
                {
                    long j = k + step;
                    sb.Append(Sub(text, k, j));
                    k = j;
                }
                else
                {
                    sb.Append('-', (int)gaps[c]);
                }
            }

            string line = sb.ToString();
            var picked = new StringBuilder();
            foreach (int col in columns)
            {
                int idx = col < 0 ? col + line.Length : col;
                if (idx < 0 || idx >= line.Length)
                    throw new ArgumentException("string index out of range");
                picked.Append(line[idx]);
            }

            lines.Add(picked.ToString());
        }

        var (sequences, coordinates) = AlignmentCoordinateFormats.ParsePrintedAlignment(lines);
        // a record keeps its identifier and annotations, a plain sequence becomes the selected letters
        var records = selected.Select((r, k) => new BioCoordSequence
        {
            Id = alignment.Sequences[r].Id,
            Length = sequences[k].Length,
            Parts = AlignmentCoordinateFormats.SeqParts([new KeyValuePair<long, string>(0, sequences[k])], sequences[k].Length),
            Description = alignment.Sequences[r].Description,
            Annotations = alignment.Sequences[r].Annotations,
        }).ToList();
        return new BioCoordAlignment { Sequences = records, Coordinates = coordinates };
    }
}
