using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentOperations
{
    private static char StrandOf(long[][] coordinates, string which)
    {
        var steps = StepsOf(coordinates);
        int m = steps[0].Length;
        bool nonNegative = true, nonPositive = true;
        for (int k = 0; k < m; k++)
        {
            int product = Math.Sign(steps[0][k]) * Math.Sign(steps[1][k]);
            if (product < 0)
                nonNegative = false;
            if (product > 0)
                nonPositive = false;
        }

        if (nonNegative)
            return '+';
        if (nonPositive)
            return '-';
        throw Bad($"Inconsistent steps in the {which} alignment");
    }

    private static void RequirePairwise(BioCoordAlignment a, string what)
    {
        RequireCoordinates(a);
        if (a.Sequences.Count != 2)
            throw Bad($"self.{what} is defined for pairwise alignments only (found alignment of {a.Sequences.Count} sequences)");
    }

    /// <summary>
    /// <c>alignment1.map(alignment2)</c>: the query of <paramref name="first"/> is the target of <paramref name="second"/> (an alignment of the reads to a transcript and the transcript to a genome give the reads on the genome, as
    /// the SAM / PSL liftover of exons); returns the alignment of <c>first.target</c> and <c>second.query</c>. The strands of the two alignments (decreasing coordinates) are combined.
    /// </summary>
    /// <exception cref="ArgumentException">The query of the first alignment and the target of the second have different lengths, the steps are inconsistent or unequal.</exception>
    public static BioCoordAlignment Map(BioCoordAlignment first, BioCoordAlignment second)
    {
        RequirePairwise(first, "query");
        RequirePairwise(second, "target");
        if (first.Sequences[1].Length != second.Sequences[0].Length)
            throw Bad($"length of alignment1 query sequence ({first.Sequences[1].Length}) != length of alignment2 target sequence ({second.Sequences[0].Length})");
        var target = first.Sequences[0];
        var query = second.Sequences[1];
        var coordinates1 = first.Coordinates.Select(r => (long[])r.Clone()).ToArray();
        var coordinates2 = second.Coordinates.Select(r => (long[])r.Clone()).ToArray();
        long n1 = first.Sequences[1].Length, n2 = second.Sequences[1].Length;
        char strand1 = StrandOf(coordinates1, "first");
        char strand2 = StrandOf(coordinates2, "second");
        if (strand1 == '+')
        {
            if (strand2 == '-')
                for (int k = 0; k < coordinates2[1].Length; k++)
                    coordinates2[1][k] = n2 - coordinates2[1][k];
        }
        else
        {
            for (int k = 0; k < coordinates1[1].Length; k++)
                coordinates1[1][k] = n1 - coordinates1[1][k];
            var original = coordinates2.Select(r => (long[])r.Clone()).ToArray();
            int m2 = original[0].Length;
            for (int k = 0; k < m2; k++)
                coordinates2[0][k] = n1 - original[0][m2 - 1 - k];
            for (int k = 0; k < m2; k++)
                coordinates2[1][k] = strand2 == '+' ? n2 - original[1][m2 - 1 - k] : original[1][m2 - 1 - k];
        }

        try
        {
            Gaps(StepsOf(coordinates1));
        }
        catch (ArgumentException)
        {
            throw Bad("Unequal step sizes in first alignment");
        }

        try
        {
            Gaps(StepsOf(coordinates2));
        }
        catch (ArgumentException)
        {
            throw Bad("Unequal step sizes in second alignment");
        }

        var path = new List<long[]>();
        long tEnd = long.MaxValue, qEnd = long.MaxValue;
        int cursor1 = 0;
        int columns1 = coordinates1[0].Length;
        long tStart1 = long.MaxValue, qStart1 = long.MaxValue, tEnd1 = 0, qEnd1 = 0;
        bool Next1()
        {
            if (cursor1 >= columns1)
                return false;
            tEnd1 = coordinates1[0][cursor1];
            qEnd1 = coordinates1[1][cursor1];
            cursor1++;
            return true;
        }

        while (Next1())
        {
            if (tStart1 < tEnd1 && qStart1 < qEnd1)
                break;
            tStart1 = tEnd1;
            qStart1 = qEnd1;
        }

        long tStart2 = long.MaxValue, qStart2 = long.MaxValue;
        for (int c2 = 0; c2 < coordinates2[0].Length; c2++)
        {
            long tEnd2 = coordinates2[0][c2], qEnd2 = coordinates2[1][c2];
            while (qStart2 < qEnd2 && tStart2 < tEnd2)
            {
                long size;
                while (true)
                {
                    if (tStart2 < qStart1)
                    {
                        size = tEnd2 < qStart1 ? tEnd2 - tStart2 : qStart1 - tStart2;
                        break;
                    }

                    if (tStart2 < qEnd1)
                    {
                        long offset = tStart2 - qStart1;
                        size = tEnd2 < qEnd1 ? tEnd2 - tStart2 : qEnd1 - tStart2;
                        long qStart = qStart2;
                        long tStart = tStart1 + offset;
                        if (tStart != tEnd || qStart != qEnd)
                        {
                            if (tStart > tEnd && qStart > qEnd)
                                path.Add([tStart, qEnd]);
                            path.Add([tStart, qStart]);
                        }

                        qEnd = qStart2 + size;
                        tEnd = tStart + size;
                        path.Add([tEnd, qEnd]);
                        break;
                    }

                    tStart1 = tEnd1;
                    qStart1 = qEnd1;
                    bool broke = false;
                    while (Next1())
                    {
                        if (tStart1 < tEnd1 && qStart1 < qEnd1)
                        {
                            broke = true;
                            break;
                        }

                        tStart1 = tEnd1;
                        qStart1 = qEnd1;
                    }

                    if (!broke)
                    {
                        size = qEnd2 - qStart2;
                        break;
                    }
                }

                qStart2 += size;
                tStart2 += size;
            }

            tStart2 = tEnd2;
            qStart2 = qEnd2;
        }

        long[][] coordinates = path.Count == 0 ? [[], []] : [[.. path.Select(p => p[0])], [.. path.Select(p => p[1])]];
        if (strand1 != strand2)
        {
            if (path.Count == 0)
                throw Bad("The mapped alignment is empty.");
            for (int k = 0; k < coordinates[1].Length; k++)
                coordinates[1][k] = n2 - coordinates[1][k];
        }

        return new BioCoordAlignment { Sequences = [target, query], Coordinates = coordinates };
    }

    /// <summary>
    /// <c>alignment.mapall(alignments)</c>: <paramref name="alignment"/> is the multiple alignment of a set of sequences (its coordinate rows are the starts of the sequences in the columns), <paramref name="alignments"/> hold for each of its
    /// rows the alignment of that sequence to a genome (or protein) sequence; the result is the multiple alignment of those genome sequences (every alignment is mapped through the multiple alignment, a factor of 3 turns protein columns into
    /// nucleotide columns).
    /// </summary>
    public static BioCoordAlignment MapAll(BioCoordAlignment alignment, IReadOnlyList<BioCoordAlignment> alignments)
    {
        RequireCoordinates(alignment);
        ArgumentNullException.ThrowIfNull(alignments);
        var list = alignments.ToList();
        long? factor = null;
        foreach (var a in list)
        {
            RequireCoordinates(a);
            var steps = StepsOf(a.Coordinates);
            var aligned = AlignedBlocks(steps);
            long step1 = 0, step2 = 0;
            for (int k = 0; k < aligned.Length; k++)
            {
                if (!aligned[k])
                    continue;
                step1 += steps[0][k];
                step2 += steps[1][k];
            }

            long step;
            if (step1 == step2 || step1 == -step2)
                step = 1;
            else if (step2 == 3 * step1)
                step = 3;
            else
                throw Bad($"unexpected steps {step1}, {step2}");
            if (factor is null)
                factor = step;
            else if (factor != step)
                throw Bad("inconsistent step sizes in alignments");
        }

        int columns = alignment.Coordinates[0].Length;
        var steps0 = new long[columns - 1];
        for (int k = 0; k < steps0.Length; k++)
        {
            long max = long.MinValue;
            foreach (var row in alignment.Coordinates)
                max = Math.Max(max, Math.Abs(row[k + 1] - row[k]));
            steps0[k] = Math.Max(0, max);
        }

        long f = factor ?? 1;
        var coordinates = new[] { new long[columns], new long[columns] };
        long cumulative = 0;
        for (int k = 1; k < columns; k++)
        {
            cumulative += steps0[k - 1];
            coordinates[0][k] = f * cumulative;
        }

        var sequence0 = new BioCoordSequence { Length = coordinates[0][^1], Parts = AlignmentCoordinateFormats.UndefinedParts(coordinates[0][^1]) };
        var mapped = new List<BioCoordAlignment>();
        for (int i = 0; i < list.Count; i++)
        {
            for (int k = 0; k < columns; k++)
                coordinates[1][k] = f * alignment.Coordinates[i][k];
            var sequence1 = new BioCoordSequence { Length = coordinates[1][^1], Parts = AlignmentCoordinateFormats.UndefinedParts(coordinates[1][^1]) };
            var alignment1 = new BioCoordAlignment { Sequences = [sequence0, sequence1], Coordinates = [(long[])coordinates[0].Clone(), (long[])coordinates[1].Clone()] };
            var coordinates2 = list[i].Coordinates.Select(r => (long[])r.Clone()).ToArray();
            for (int k = 0; k < coordinates2[0].Length; k++)
                coordinates2[0][k] *= f;
            var alignment2 = new BioCoordAlignment { Sequences = [sequence1, list[i].Sequences[1]], Coordinates = coordinates2 };
            mapped.Add(Map(alignment1, alignment2));
        }

        var rowsOut = mapped.Select(_ => new List<long>()).ToArray();
        var remaining = mapped.Select(a => a.Coordinates.Select(r => r.ToList()).ToArray()).ToArray();
        bool done = false;
        long previous = 0;
        while (!done)
        {
            done = true;
            long position = remaining.Where(c => c[0].Count > 0).Min(c => c[0][0]);
            for (int i = 0; i < remaining.Length; i++)
            {
                var c = remaining[i];
                if (c[0].Count == 0)
                {
                    rowsOut[i].Add(rowsOut[i][^1]);
                }
                else if (c[0][0] == position)
                {
                    rowsOut[i].Add(c[1][0]);
                    c[0].RemoveAt(0);
                    c[1].RemoveAt(0);
                    if (c[0].Any(x => x != 0) || c[1].Any(x => x != 0))
                        done = false;
                }
                else if (c[0][0] > position)
                {
                    if (rowsOut[i].Count > 0)
                    {
                        long step = c[1][0] > rowsOut[i][^1] ? position - previous : 0;
                        rowsOut[i].Add(rowsOut[i][^1] + step);
                    }
                    else
                    {
                        rowsOut[i].Add(c[1][0]);
                    }
                }
                else
                {
                    throw new InvalidOperationException("mapall: inconsistent positions.");
                }
            }

            previous = position;
        }

        return new BioCoordAlignment { Sequences = [.. mapped.Select(a => a.Sequences[1])], Coordinates = [.. rowsOut.Select(r => r.ToArray())] };
    }
}
