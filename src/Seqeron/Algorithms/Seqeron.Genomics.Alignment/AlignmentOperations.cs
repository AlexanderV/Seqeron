using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>A Python slice <c>start:stop:step</c> (null = omitted) for the index arguments of <see cref="AlignmentOperations"/> (<c>alignment[rows, cols]</c>).</summary>
/// <param name="Start">The start index (negative counts from the end) or null.</param>
/// <param name="Stop">The stop index (exclusive) or null.</param>
/// <param name="Step">The step (default 1) or null.</param>
public readonly record struct PySlice(int? Start = null, int? Stop = null, int? Step = null)
{
    /// <summary>The slice <c>:</c> (everything).</summary>
    public static PySlice All => new();

    /// <summary>Python's <c>slice.indices(length)</c>.</summary>
    public (int Start, int Stop, int Step) Indices(int length)
    {
        int step = Step ?? 1;
        if (step == 0)
            throw new ArgumentException("slice step cannot be zero");
        int start, stop;
        if (step > 0)
        {
            start = Start is null ? 0 : Clamp(Start.Value, 0, length);
            stop = Stop is null ? length : Clamp(Stop.Value, 0, length);
        }
        else
        {
            start = Start is null ? length - 1 : Clamp(Start.Value, -1, length - 1);
            stop = Stop is null ? -1 : Clamp(Stop.Value, -1, length - 1);
        }

        return (start, stop, step);

        int Clamp(int v, int lo, int hi)
        {
            if (v < 0)
            {
                v += length;
                if (v < lo)
                    v = lo;
            }
            else if (v > hi)
            {
                v = hi;
            }

            return v;
        }
    }

    /// <summary>The indices the slice selects of a sequence of <paramref name="length"/> items.</summary>
    public IReadOnlyList<int> Range(int length)
    {
        var (start, stop, step) = Indices(length);
        var list = new List<int>();
        if (step > 0)
        {
            for (int i = start; i < stop; i += step)
                list.Add(i);
        }
        else
        {
            for (int i = start; i > stop; i += step)
                list.Add(i);
        }

        return list;
    }
}

/// <summary>The substitution counts of <see cref="AlignmentOperations.Substitutions"/> (Biopython's <c>alignment.substitutions</c>, a 2-D <c>substitution_matrices.Array</c>): <c>Values[i, j]</c> counts the pairs (letter i of an earlier sequence, letter j of a later one) aligned to each other.</summary>
/// <param name="Letters">The sorted letters of all sequences (rows and columns).</param>
/// <param name="Values">The counts.</param>
public sealed record SubstitutionCounts(string Letters, double[,] Values)
{
    /// <summary>The count of the pair (<paramref name="a"/>, <paramref name="b"/>).</summary>
    public double this[char a, char b] => Values[Letters.IndexOf(a, StringComparison.Ordinal), Letters.IndexOf(b, StringComparison.Ordinal)];
}

/// <summary>
/// The operations of Biopython 1.88's <c>Bio.Align.Alignment</c> that work on the block coordinates of an alignment of any number of sequences (<see cref="BioCoordAlignment"/>, the object the coordinate-format readers
/// and <see cref="AlignmentResultCoordinates.ToBioCoordAlignment(AlignmentResult, string, string, string, string, AlignmentBridgeOptions?)"/> produce): shape / length, rows and columns (<c>alignment[i]</c>, <c>alignment[i, j]</c>, slices,
/// column iterables), <c>aligned</c>, <c>indices</c>, <c>inverse_indices</c>, <c>frequencies</c>, <c>substitutions</c>, <c>counts</c> (<see cref="AlignmentCounts"/>), <c>reverse_complement</c>, <c>sort</c>, <c>+</c>, the ordering of alignments,
/// <c>map</c>, <c>mapall</c>, <c>from_alignments_with_same_reference</c>, <c>parse_printed_alignment</c> and the pretty text <c>str(alignment)</c>. Every method is a port of the Python code (same errors for the same inputs); where Biopython changes the
/// object in place (<c>sort</c>) a new alignment is returned.
/// </summary>
public static partial class AlignmentOperations
{
    private static ArgumentException Bad(string message) => new(message);

    private sealed record Oriented(BioCoordSequence[] Sequences, long[][] Coordinates, long[][] Steps);

    private static long[][] StepsOf(long[][] c)
    {
        var steps = new long[c.Length][];
        for (int i = 0; i < c.Length; i++)
        {
            steps[i] = new long[Math.Max(0, c[i].Length - 1)];
            for (int k = 0; k < steps[i].Length; k++)
                steps[i][k] = c[i][k + 1] - c[i][k];
        }

        return steps;
    }

    // sum(steps != 0, 0) > 1: the blocks in which at least two sequences advance
    private static bool[] AlignedBlocks(long[][] steps)
    {
        int m = steps.Length == 0 ? 0 : steps[0].Length;
        var aligned = new bool[m];
        for (int k = 0; k < m; k++)
        {
            int count = 0;
            foreach (var row in steps)
                if (row[k] != 0)
                    count++;
            aligned[k] = count > 1;
        }

        return aligned;
    }

    // the strand of a row from its aligned steps: +1 all >= 0, -1 all <= 0 (and not all >= 0), else inconsistent
    private static int Strand(long[] row, bool[] aligned, int index)
    {
        bool nonNegative = true, nonPositive = true;
        for (int k = 0; k < row.Length; k++)
        {
            if (!aligned[k])
                continue;
            if (row[k] < 0)
                nonNegative = false;
            if (row[k] > 0)
                nonPositive = false;
        }

        if (nonNegative)
            return 1;
        if (nonPositive)
            return -1;
        throw Bad($"Inconsistent steps in row {index}");
    }

    private static void RequireCoordinates(BioCoordAlignment a)
    {
        ArgumentNullException.ThrowIfNull(a);
        if (a.NoCoordinates)
            throw Bad("The alignment has no coordinates.");
    }

    // the sequences with the reverse-strand rows reverse complemented, their coordinates mirrored and the steps made non-negative (the preamble of __array__, frequencies, __getitem__, substitutions)
    private static Oriented Orient(BioCoordAlignment a)
    {
        RequireCoordinates(a);
        var sequences = a.Sequences.ToArray();
        var coordinates = a.Coordinates.Select(r => (long[])r.Clone()).ToArray();
        var steps = StepsOf(coordinates);
        var aligned = AlignedBlocks(steps);
        for (int i = 0; i < sequences.Length; i++)
        {
            if (Strand(steps[i], aligned, i) < 0)
            {
                long length = sequences[i].Length;
                sequences[i] = AlignmentCoordinateFormats.ReverseComplement(sequences[i]);
                for (int k = 0; k < coordinates[i].Length; k++)
                    coordinates[i][k] = length - coordinates[i][k];
                for (int k = 0; k < steps[i].Length; k++)
                    steps[i][k] = -steps[i][k];
            }
        }

        return new Oriented(sequences, coordinates, steps);
    }

    private static long[] Gaps(long[][] steps)
    {
        int m = steps.Length == 0 ? 0 : steps[0].Length;
        var gaps = new long[m];
        for (int k = 0; k < m; k++)
        {
            long max = long.MinValue;
            foreach (var row in steps)
                max = Math.Max(max, row[k]);
            gaps[k] = max;
        }

        foreach (var row in steps)
            for (int k = 0; k < m; k++)
                if (row[k] != gaps[k] && row[k] > 0)
                    throw Bad("Unequal step sizes in alignment");
        return gaps;
    }

    private static string Contents(BioCoordSequence s) => s.Length < 0 ? throw Bad("The sequence is undefined (UndefinedSequenceError).") : s.Slice(0, s.Length) ?? throw Bad("The sequence contents are undefined (UndefinedSequenceError).");

    private static string Sub(string s, long start, long end)
    {
        long a = Math.Max(0, Math.Min(start, s.Length)), b = Math.Max(0, Math.Min(end, s.Length));
        return b <= a ? "" : s.Substring((int)a, (int)(b - a));
    }

    private static BioCoordSequence Plain(string text, string? id = null) => new() { Id = id, Length = text.Length, Parts = AlignmentCoordinateFormats.SeqParts([new KeyValuePair<long, string>(0, text)], text.Length) };

    /// <summary><c>len(alignment)</c>: the number of sequences.</summary>
    public static int Count(BioCoordAlignment alignment)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        return alignment.Sequences.Count;
    }

    /// <summary><c>alignment.length</c>: the number of columns of the printed alignment.</summary>
    public static int Length(BioCoordAlignment alignment)
    {
        RequireCoordinates(alignment);
        int n = alignment.Coordinates.Length;
        if (n == 0)
            return 0;
        var steps = StepsOf(alignment.Coordinates);
        var aligned = AlignedBlocks(steps);
        for (int i = 0; i < n; i++)
        {
            if (Strand(steps[i], aligned, i) < 0)
                for (int k = 0; k < steps[i].Length; k++)
                    steps[i][k] = -steps[i][k];
        }

        return (int)Gaps(steps).Sum();
    }

    /// <summary><c>alignment.shape</c>: (number of rows, number of columns).</summary>
    public static (int Rows, int Columns) Shape(BioCoordAlignment alignment) => (RequireCoordinatesCount(alignment), Length(alignment));

    private static int RequireCoordinatesCount(BioCoordAlignment a)
    {
        RequireCoordinates(a);
        return a.Coordinates.Length;
    }

    /// <summary><c>alignment[index]</c>: the printed row of a sequence (negative = from the end).</summary>
    public static string Row(BioCoordAlignment alignment, int index)
    {
        RequireCoordinates(alignment);
        int n = alignment.Coordinates.Length;
        if (index < 0)
            index += n;
        if (index < 0 || index >= n)
            throw new ArgumentException("row index out of range");
        return AlignmentCoordinateFormats.GetRow(alignment, index);
    }

    /// <summary>The printed rows of every sequence (what <c>numpy.array(alignment)</c> holds, row by row; <c>__array__</c>).</summary>
    public static string[] Rows(BioCoordAlignment alignment)
    {
        var o = Orient(alignment);
        var gaps = Gaps(o.Steps);
        var rows = new string[o.Sequences.Length];
        for (int i = 0; i < rows.Length; i++)
        {
            var sb = new StringBuilder();
            string text = Contents(o.Sequences[i]);
            long k = o.Coordinates[i][0];
            for (int c = 0; c < gaps.Length; c++)
            {
                long step = o.Steps[i][c];
                if (step > 0)
                {
                    sb.Append(Sub(text, k, k + step));
                    k += step;
                }
                else if (step < 0)
                {
                    k += step;
                }
                else
                {
                    sb.Append('-', (int)gaps[c]);
                }
            }

            rows[i] = sb.ToString();
        }

        return rows;
    }

    /// <summary><c>alignment.aligned</c>: for each of the two sequences the (start, end) of every block in which both are aligned (a reverse-strand row has start &gt; end).</summary>
    /// <returns><c>result[row][block] = [start, end]</c>.</returns>
    public static long[][][] Aligned(BioCoordAlignment alignment)
    {
        RequireCoordinates(alignment);
        if (alignment.Sequences.Count > 2)
            throw new NotSupportedException("aligned is currently implemented for pairwise alignments only");
        var coordinates = alignment.Coordinates.Select(r => (long[])r.Clone()).ToArray();
        var steps = StepsOf(coordinates);
        var aligned = AlignedBlocks(steps);
        var mirrored = new bool[coordinates.Length];
        for (int i = 0; i < coordinates.Length; i++)
        {
            if (Strand(steps[i], aligned, i) < 0)
            {
                mirrored[i] = true;
                long len = alignment.Sequences[i].Length;
                for (int k = 0; k < coordinates[i].Length; k++)
                    coordinates[i][k] = len - coordinates[i][k];
                for (int k = 0; k < steps[i].Length; k++)
                    steps[i][k] = -steps[i][k];
            }
        }

        int m = coordinates.Length == 0 ? 0 : coordinates[0].Length;
        var segments = coordinates.Select(_ => new List<long[]>()).ToArray();
        for (int k = 0; k + 1 < m; k++)
        {
            long min = long.MaxValue;
            for (int i = 0; i < coordinates.Length; i++)
                min = Math.Min(min, Math.Abs(coordinates[i][k + 1] - coordinates[i][k]));
            if (min == 0)
                continue;
            for (int i = 0; i < coordinates.Length; i++)
                segments[i].Add([coordinates[i][k], coordinates[i][k + 1]]);
        }

        var result = new long[coordinates.Length][][];
        for (int i = 0; i < result.Length; i++)
        {
            long len = alignment.Sequences[i].Length;
            result[i] = [.. segments[i].Select(s => mirrored[i] ? new[] { len - s[0], len - s[1] } : s)];
        }

        return result;
    }

    /// <summary><c>alignment.indices</c>: for every column the index in each sequence of the letter in it (-1 for a gap); a reverse-strand row counts down.</summary>
    /// <returns><c>result[row][column]</c>.</returns>
    public static int[][] Indices(BioCoordAlignment alignment)
    {
        RequireCoordinates(alignment);
        var (n, m) = Shape(alignment);
        var a = Enumerable.Range(0, n).Select(_ => Enumerable.Repeat(-1, m).ToArray()).ToArray();
        var coordinates = alignment.Coordinates;
        var steps = StepsOf(coordinates);
        var aligned = AlignedBlocks(steps);
        var rcs = new bool[n];
        for (int i = 0; i < n; i++)
            rcs[i] = Strand(steps[i], aligned, i) < 0;
        int col = 0, next = 0;
        for (int k = 1; k < coordinates[0].Length; k++)
        {
            for (int r = 0; r < n; r++)
            {
                long start = coordinates[r][k - 1], end = coordinates[r][k];
                if (!rcs[r] && start < end)
                {
                    next = (int)(col + end - start);
                    for (int t = 0; t < end - start; t++)
                        a[r][col + t] = (int)(start + t);
                }
                else if (rcs[r] && start > end)
                {
                    next = (int)(col + start - end);
                    for (int t = 0; t < start - end; t++)
                        a[r][col + t] = (int)(start - 1 - t);
                }
            }

            col = next;
        }

        return a;
    }

    /// <summary><c>alignment.inverse_indices</c>: for every letter of every sequence the column it is in (-1 for a letter that is not aligned).</summary>
    /// <returns><c>result[row][letter]</c>.</returns>
    public static int[][] InverseIndices(BioCoordAlignment alignment)
    {
        RequireCoordinates(alignment);
        int n = alignment.Coordinates.Length;
        var a = alignment.Sequences.Select(s => Enumerable.Repeat(-1, (int)Math.Max(0, s.Length)).ToArray()).ToArray();
        var coordinates = alignment.Coordinates;
        var steps = StepsOf(coordinates);
        var aligned = AlignedBlocks(steps);
        var rcs = new bool[n];
        for (int i = 0; i < n; i++)
            rcs[i] = Strand(steps[i], aligned, i) < 0;
        int col = 0, next = 0;
        for (int k = 0; k + 1 < coordinates[0].Length; k++)
        {
            for (int r = 0; r < n; r++)
            {
                long start = coordinates[r][k], end = coordinates[r][k + 1];
                if (!rcs[r] && start < end)
                {
                    next = (int)(col + end - start);
                    for (int t = 0; t < end - start; t++)
                        a[r][start + t] = col + t;
                }
                else if (rcs[r] && start > end)
                {
                    next = (int)(col + start - end);
                    for (int t = 0; t < start - end; t++)
                        a[r][start - 1 - t] = col + t;
                }
            }

            col = next;
        }

        return a;
    }

    /// <summary>
    /// <c>alignment.frequencies</c>: for every letter (and <c>-</c> for gaps) the weight of the sequences that have it in each column (the weight of a sequence is the <c>weight</c> annotation of its record, default 1).
    /// The dictionary is in order of first appearance.
    /// </summary>
    public static IReadOnlyDictionary<char, double[]> Frequencies(BioCoordAlignment alignment)
    {
        var o = Orient(alignment);
        var gaps = Gaps(o.Steps);
        int length = (int)gaps.Sum();
        var counts = new Dictionary<char, double[]>();
        double[] Row(char c)
        {
            if (!counts.TryGetValue(c, out var row))
                counts[c] = row = new double[length];
            return row;
        }

        for (int i = 0; i < o.Sequences.Length; i++)
        {
            double weight = 1.0;
            if (o.Sequences[i].Annotations.TryGetPropertyValue("weight", out var w) && w is not null)
                weight = w.GetValue<double>();
            string text = Contents(o.Sequences[i]);
            long k = o.Coordinates[i][0];
            int m = 0;
            for (int c = 0; c < gaps.Length; c++)
            {
                long step = o.Steps[i][c];
                if (step > 0)
                {
                    long j = k + step;
                    string sub = Sub(text, k, j);
                    for (int t = 0; t < sub.Length; t++)
                        Row(sub[t])[m + t] += weight;
                    m += (int)step;
                    k = j;
                }
                else if (step < 0)
                {
                    k += step;
                }
                else
                {
                    var row = Row('-');
                    for (int t = 0; t < gaps[c]; t++)
                        row[m + t] += weight;
                    m += (int)gaps[c];
                }
            }
        }

        return counts;
    }

    /// <summary><c>alignment.substitutions</c>: how often each letter is aligned to each other letter (pairs of sequences i &lt; j; the first sequence of the pair gives the row).</summary>
    public static SubstitutionCounts Substitutions(BioCoordAlignment alignment)
    {
        RequireCoordinates(alignment);
        var coordinates = alignment.Coordinates.Select(r => (long[])r.Clone()).ToArray();
        var sequences = alignment.Sequences.ToArray();
        var steps = StepsOf(alignment.Coordinates);
        var aligned = AlignedBlocks(steps);
        for (int i = 0; i < sequences.Length; i++)
        {
            if (Strand(steps[i], aligned, i) < 0)
            {
                long len = sequences[i].Length;
                sequences[i] = AlignmentCoordinateFormats.ReverseComplement(sequences[i]);
                for (int k = 0; k < coordinates[i].Length; k++)
                    coordinates[i][k] = len - coordinates[i][k];
            }
        }

        var letters = new SortedSet<char>();
        var texts = new string[sequences.Length];
        for (int i = 0; i < sequences.Length; i++)
        {
            var s = sequences[i];
            if (s.Parts is null)
                throw Bad("The sequence contents are undefined (UndefinedSequenceError).");
            foreach (var (_, part) in s.Parts)
                foreach (char c in part)
                    letters.Add(c);
            texts[i] = s.Slice(0, s.Length) ?? "";
        }

        string alphabet = new([.. letters]);
        var m = new double[alphabet.Length, alphabet.Length];
        int n = sequences.Length;
        for (int i1 = 0; i1 < n; i1++)
        {
            for (int i2 = i1 + 1; i2 < n; i2++)
            {
                long start1 = long.MaxValue, start2 = long.MaxValue;
                for (int k = 0; k < coordinates[i1].Length; k++)
                {
                    long end1 = coordinates[i1][k], end2 = coordinates[i2][k];
                    if (start1 < end1 && start2 < end2)
                    {
                        string seg1 = Sub(texts[i1], start1, end1), seg2 = Sub(texts[i2], start2, end2);
                        if (seg1.Length != seg2.Length)
                            throw Bad("Unequal step sizes in alignment");
                        for (int t = 0; t < seg1.Length; t++)
                            m[alphabet.IndexOf(seg1[t], StringComparison.Ordinal), alphabet.IndexOf(seg2[t], StringComparison.Ordinal)] += 1.0;
                    }

                    start1 = end1;
                    start2 = end2;
                }
            }
        }

        return new SubstitutionCounts(alphabet, m);
    }

    /// <summary><c>alignment.reverse_complement()</c>: every sequence reverse complemented, the coordinates mirrored and their order reversed (the score and the other attributes are dropped, as Biopython builds a new <c>Alignment</c>).</summary>
    public static BioCoordAlignment ReverseComplement(BioCoordAlignment alignment)
    {
        RequireCoordinates(alignment);
        var sequences = alignment.Sequences.Select(AlignmentCoordinateFormats.ReverseComplement).ToList();
        var coordinates = new long[sequences.Count][];
        for (int i = 0; i < coordinates.Length; i++)
            coordinates[i] = [.. alignment.Coordinates[i].Reverse().Select(x => sequences[i].Length - x)];
        return new BioCoordAlignment { Sequences = sequences, Coordinates = coordinates };
    }

    /// <summary>
    /// <c>alignment.sort(key, reverse)</c>: the sequences ordered by their identifier (or, when a sequence has none, by their residues) or by <paramref name="key"/> (stable; <paramref name="reverse"/> keeps the order of equal keys, as Python's <c>sorted</c>);
    /// the coordinate rows move with them. Biopython sorts in place; this returns the sorted alignment.
    /// </summary>
    public static BioCoordAlignment Sort(BioCoordAlignment alignment, Func<BioCoordSequence, string>? key = null, bool reverse = false)
    {
        RequireCoordinates(alignment);
        var sequences = alignment.Sequences;
        Func<BioCoordSequence, string> selector = key ?? (sequences.All(s => s.Id is not null) ? s => s.Id! : Contents);
        var keys = sequences.Select(selector).ToList();
        var order = Enumerable.Range(0, sequences.Count);
        var indices = (reverse ? order.OrderByDescending(i => keys[i], StringComparer.Ordinal) : order.OrderBy(i => keys[i], StringComparer.Ordinal)).ToList();
        return new BioCoordAlignment
        {
            Sequences = [.. indices.Select(i => sequences[i])],
            Coordinates = [.. indices.Select(i => (long[])alignment.Coordinates[i].Clone())],
            Attributes = alignment.Attributes,
        };
    }

    /// <summary>
    /// <c>alignment1 + alignment2</c>: the two alignments of the same sequences end to end (the first sequence's aligned part followed by the second's); the result keeps the annotations the two have in common
    /// (<c>attrs.annotations</c>). The joining column is kept twice in the coordinates, as Biopython does.
    /// </summary>
    public static BioCoordAlignment Add(BioCoordAlignment left, BioCoordAlignment right)
    {
        RequireCoordinates(left);
        RequireCoordinates(right);
        ArgumentNullException.ThrowIfNull(right);
        if (left.Sequences.Count != right.Sequences.Count)
            throw Bad("When adding two alignments they must have the same length (i.e. same number of rows)");
        int n = left.Sequences.Count;
        var sequences = new List<BioCoordSequence>();
        var coordinates = new long[n][];
        for (int i = 0; i < n; i++)
        {
            long start1 = left.Coordinates[i][0], end1 = left.Coordinates[i][^1], start2 = right.Coordinates[i][0], end2 = right.Coordinates[i][^1];
            sequences.Add(Plain(Sub(Contents(left.Sequences[i]), start1, end1) + Sub(Contents(right.Sequences[i]), start2, end2)));
            long offset = start2 - end1 + start1;
            coordinates[i] = [.. left.Coordinates[i].Select(x => x - start1), .. right.Coordinates[i].Select(x => x - offset)];
        }

        var attributes = new System.Text.Json.Nodes.JsonObject();
        if (left.Attributes["annotations"] is System.Text.Json.Nodes.JsonObject la && right.Attributes["annotations"] is System.Text.Json.Nodes.JsonObject ra)
        {
            var common = new System.Text.Json.Nodes.JsonObject();
            foreach (var (k, v) in la)
                if (ra.TryGetPropertyValue(k, out var other) && other?.ToJsonString() == v?.ToJsonString())
                    common[k] = v?.DeepClone();
            attributes["annotations"] = common;
        }

        return new BioCoordAlignment { Sequences = sequences, Coordinates = coordinates, Attributes = attributes };
    }

    private static int CompareSequences(BioCoordAlignment a, BioCoordAlignment b)
    {
        int n = Math.Max(a.Sequences.Count, b.Sequences.Count);
        for (int i = 0; i < n; i++)
        {
            if (i >= a.Sequences.Count || i >= b.Sequences.Count)
                throw Bad("'<' not supported between an alignment sequence and None");
            int c = string.CompareOrdinal(Contents(a.Sequences[i]), Contents(b.Sequences[i]));
            if (c != 0)
                return c < 0 ? -1 : 1;
        }

        return 0;
    }

    /// <summary>
    /// The order of <c>alignment1 &lt; alignment2</c> (<c>__lt__</c> and friends): the sequences compared as strings one after the other, then the coordinate columns as tuples. Returns -1, 0 or 1;
    /// alignments with the same sequences and coordinates are equal. Used by <c>sorted(alignments)</c>.
    /// </summary>
    public static int Compare(BioCoordAlignment left, BioCoordAlignment right)
    {
        RequireCoordinates(left);
        RequireCoordinates(right);
        int c = CompareSequences(left, right);
        if (c != 0)
            return c;
        int columns = Math.Min(left.Coordinates.Length == 0 ? 0 : left.Coordinates[0].Length, right.Coordinates.Length == 0 ? 0 : right.Coordinates[0].Length);
        int rows = Math.Min(left.Coordinates.Length, right.Coordinates.Length);
        for (int k = 0; k < columns; k++)
        {
            for (int r = 0; r < rows; r++)
            {
                long l = left.Coordinates[r][k], rr = right.Coordinates[r][k];
                if (l != rr)
                    return l < rr ? -1 : 1;
            }

            if (left.Coordinates.Length != right.Coordinates.Length)
                return left.Coordinates.Length < right.Coordinates.Length ? -1 : 1;
        }

        return 0;
    }

    /// <summary><c>alignment1 == alignment2</c>: the same sequences and the same coordinates (<c>numpy.array_equal</c>).</summary>
    public static bool AreEqual(BioCoordAlignment left, BioCoordAlignment right)
    {
        RequireCoordinates(left);
        RequireCoordinates(right);
        int n = Math.Max(left.Sequences.Count, right.Sequences.Count);
        for (int i = 0; i < n; i++)
        {
            if (i >= left.Sequences.Count || i >= right.Sequences.Count)
                return false;
            if (!string.Equals(Contents(left.Sequences[i]), Contents(right.Sequences[i]), StringComparison.Ordinal))
                return false;
        }

        if (left.Coordinates.Length != right.Coordinates.Length)
            return false;
        for (int i = 0; i < left.Coordinates.Length; i++)
            if (!left.Coordinates[i].SequenceEqual(right.Coordinates[i]))
                return false;
        return true;
    }

    /// <summary><c>Alignment.parse_printed_alignment(lines)</c>: the ungapped sequences and the block coordinates of equally wide printed rows (gaps are <c>-</c>).</summary>
    public static (IReadOnlyList<string> Sequences, long[][] Coordinates) ParsePrintedAlignment(IReadOnlyList<string> lines)
    {
        ArgumentNullException.ThrowIfNull(lines);
        var (sequences, coordinates) = AlignmentCoordinateFormats.ParsePrintedAlignment(lines);
        return (sequences, coordinates);
    }

    /// <summary>
    /// <c>Alignment.from_alignments_with_same_reference(alignments)</c>: pairwise (or multiple) alignments whose first sequence is the same reference, combined into one multiple alignment (the reference, then the other
    /// sequences of every alignment in order).
    /// </summary>
    /// <exception cref="ArgumentException">No alignment, different lengths / contents of the references, inconsistent reference coordinates.</exception>
    public static BioCoordAlignment FromAlignmentsWithSameReference(IReadOnlyList<BioCoordAlignment> alignments)
    {
        ArgumentNullException.ThrowIfNull(alignments);
        if (alignments.Count == 0)
            throw Bad("No pairwise alignments provided.");
        foreach (var alignment in alignments)
            RequireCoordinates(alignment);
        var firstSeqs = alignments.Select(a => a.Sequences[0]).ToList();
        if (!firstSeqs.All(s => s.Length == firstSeqs[0].Length))
            throw Bad("All reference sequences must have the same length.");
        var defined = new HashSet<string>(StringComparer.Ordinal);
        foreach (var s in firstSeqs)
        {
            string? text = s.Length < 0 ? null : s.Slice(0, s.Length);
            if (text is not null)
                defined.Add(text.ToUpperInvariant());
        }

        if (defined.Count > 1)
            throw Bad("All reference sequences must match (excluding gaps).");
        var sequences = new List<BioCoordSequence> { alignments[0].Sequences[0] };
        sequences.AddRange(alignments.SelectMany(a => a.Sequences.Skip(1)));
        // the coordinates of reference vs every other sequence, one (m x 2) array each
        var paired = new List<(long[] Reference, long[] Query)>();
        foreach (var alignment in alignments)
            for (int row = 1; row < alignment.Coordinates.Length; row++)
                paired.Add((alignment.Coordinates[0], alignment.Coordinates[row]));
        if (paired.Select(p => p.Reference[0]).Distinct().Count() != 1 || paired.Select(p => p.Reference[^1]).Distinct().Count() != 1)
            throw Bad("Reference coordinates do not align consistently across pairwise alignments.");
        long referencePosition = paired[0].Reference[0];
        int q = paired.Count;
        var cursor = new int[q];
        var positions = new long[q][];
        var nextPositions = new long[q][];
        for (int i = 0; i < q; i++)
            positions[i] = [paired[i].Reference[cursor[i]], paired[i].Query[cursor[i]++]];
        for (int i = 0; i < q; i++)
            nextPositions[i] = Next(i);
        var msa = new List<long[]>();
        msa.Add([referencePosition, .. positions.Select(p => p[1])]);
        while (true)
        {
            if (nextPositions.All(p => p[0] == long.MaxValue))
                break;
            var targetSteps = nextPositions.Select(p => p[0] == long.MaxValue ? long.MaxValue - referencePosition : p[0] - referencePosition).ToArray();
            int index = 0;
            for (int i = 1; i < q; i++)
                if (targetSteps[i] < targetSteps[index])
                    index = i;
            long targetStep = targetSteps[index];
            var queryStepsAll = new long[q];
            for (int i = 0; i < q; i++)
                queryStepsAll[i] = nextPositions[i][1] - positions[i][1];
            long queryStep = queryStepsAll[index];
            referencePosition += targetStep;
            positions[index] = [nextPositions[index][0], nextPositions[index][1]];
            nextPositions[index] = Next(index);
            if (targetStep == 0)
            {
                if (queryStep == 0)
                    continue;
            }
            else
            {
                for (int i = 0; i < q; i++)
                {
                    if (i == index)
                        continue;
                    if (queryStepsAll[i] > 0)
                    {
                        positions[i][0] += targetStep;
                        positions[i][1] += targetStep;
                    }
                    else
                    {
                        positions[i][0] += targetStep;
                    }
                }
            }

            msa.Add([referencePosition, .. positions.Select(p => p[1])]);
        }

        int rowsOut = q + 1;
        var coordinates = new long[rowsOut][];
        for (int r = 0; r < rowsOut; r++)
            coordinates[r] = [.. msa.Select(col => col[r])];
        return new BioCoordAlignment { Sequences = sequences, Coordinates = coordinates };

        long[] Next(int i) => cursor[i] < paired[i].Reference.Length ? [paired[i].Reference[cursor[i]], paired[i].Query[cursor[i]++]] : [long.MaxValue, long.MaxValue];
    }

    /// <summary>The text of <see cref="Format(BioCoordAlignment)"/> as <c>str(alignment)</c>.</summary>
    public static string ToText(BioCoordAlignment alignment) => Format(alignment);

    internal static string FormatInvariant(long v) => v.ToString(CultureInfo.InvariantCulture);
}
