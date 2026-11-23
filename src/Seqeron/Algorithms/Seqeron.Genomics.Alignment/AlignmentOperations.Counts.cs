using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The counts of an alignment (Biopython 1.88 <c>AlignmentCounts</c>, <c>alignment.counts(scoring)</c>): aligned letters, identities, mismatches, positives and the gaps split by side (left / internal / right), direction
/// (insertion = a gap in the target, deletion = a gap in the query) and opened / extended, with the substitution score, the gap score and the total score when a scoring is given. The counts add up over every pair of sequences of the alignment.
/// </summary>
public sealed record AlignmentCounts
{
    /// <summary>Insertion gaps opened on the left of the alignment.</summary>
    public long OpenLeftInsertions { get; init; }

    /// <summary>Insertion gap extensions on the left.</summary>
    public long ExtendLeftInsertions { get; init; }

    /// <summary>Deletion gaps opened on the left.</summary>
    public long OpenLeftDeletions { get; init; }

    /// <summary>Deletion gap extensions on the left.</summary>
    public long ExtendLeftDeletions { get; init; }

    /// <summary>Insertion gaps opened inside the alignment.</summary>
    public long OpenInternalInsertions { get; init; }

    /// <summary>Insertion gap extensions inside.</summary>
    public long ExtendInternalInsertions { get; init; }

    /// <summary>Deletion gaps opened inside.</summary>
    public long OpenInternalDeletions { get; init; }

    /// <summary>Deletion gap extensions inside.</summary>
    public long ExtendInternalDeletions { get; init; }

    /// <summary>Insertion gaps opened on the right.</summary>
    public long OpenRightInsertions { get; init; }

    /// <summary>Insertion gap extensions on the right.</summary>
    public long ExtendRightInsertions { get; init; }

    /// <summary>Deletion gaps opened on the right.</summary>
    public long OpenRightDeletions { get; init; }

    /// <summary>Deletion gap extensions on the right.</summary>
    public long ExtendRightDeletions { get; init; }

    /// <summary>Number of aligned letters (pairs).</summary>
    public long Aligned { get; init; }

    /// <summary>Number of identical aligned letters.</summary>
    public long Identities { get; init; }

    /// <summary>Number of different aligned letters.</summary>
    public long Mismatches { get; init; }

    /// <summary>Number of aligned letters with a positive substitution score; null without a substitution matrix.</summary>
    public long? Positives { get; init; }

    /// <summary>The total gap score; null when no aligner was given.</summary>
    public double? GapScore { get; init; }

    /// <summary>The total substitution score; null when neither an aligner nor a matrix was given.</summary>
    public double? SubstitutionScore { get; init; }

    /// <summary>The alignment score (gap score + substitution score); null when either is unknown.</summary>
    public double? Score => GapScore is { } g && SubstitutionScore is { } s ? g + s : null;

    /// <summary>Letters inserted on the left (opened + extended).</summary>
    public long LeftInsertions => OpenLeftInsertions + ExtendLeftInsertions;

    /// <summary>Letters deleted on the left.</summary>
    public long LeftDeletions => OpenLeftDeletions + ExtendLeftDeletions;

    /// <summary>Gap positions on the left.</summary>
    public long LeftGaps => LeftInsertions + LeftDeletions;

    /// <summary>Letters inserted inside.</summary>
    public long InternalInsertions => OpenInternalInsertions + ExtendInternalInsertions;

    /// <summary>Letters deleted inside.</summary>
    public long InternalDeletions => OpenInternalDeletions + ExtendInternalDeletions;

    /// <summary>Gap positions inside.</summary>
    public long InternalGaps => InternalInsertions + InternalDeletions;

    /// <summary>Letters inserted on the right.</summary>
    public long RightInsertions => OpenRightInsertions + ExtendRightInsertions;

    /// <summary>Letters deleted on the right.</summary>
    public long RightDeletions => OpenRightDeletions + ExtendRightDeletions;

    /// <summary>Gap positions on the right.</summary>
    public long RightGaps => RightInsertions + RightDeletions;

    /// <summary>Insertion positions in total.</summary>
    public long Insertions => LeftInsertions + InternalInsertions + RightInsertions;

    /// <summary>Deletion positions in total.</summary>
    public long Deletions => LeftDeletions + InternalDeletions + RightDeletions;

    /// <summary>Gap positions in total.</summary>
    public long Gaps => LeftGaps + InternalGaps + RightGaps;

    /// <summary>Gaps opened in total.</summary>
    public long OpenGaps => OpenLeftInsertions + OpenLeftDeletions + OpenInternalInsertions + OpenInternalDeletions + OpenRightInsertions + OpenRightDeletions;

    /// <summary>Gap extensions in total.</summary>
    public long ExtendGaps => ExtendLeftInsertions + ExtendLeftDeletions + ExtendInternalInsertions + ExtendInternalDeletions + ExtendRightInsertions + ExtendRightDeletions;

    private static string Real(double d) => AlignmentCoordinateFormats.PyFloat(d);

    /// <summary>The text <c>str(counts)</c> of Biopython prints.</summary>
    public override string ToString()
    {
        double score = (GapScore ?? double.NaN) + (SubstitutionScore ?? double.NaN);
        var sb = new StringBuilder("AlignmentCounts object with\n");
        var sub = SubstitutionScore;
        var gap = GapScore;
        if (!double.IsNaN(score))
        {
            sb.Append("    score = ").Append(Real(score)).Append(":\n");
            if (sub is not null)
                sb.Append("        substitution_score = ").Append(Real(sub.Value)).Append(",\n");
            if (gap is not null)
                sb.Append("        gap_score = ").Append(Real(gap.Value)).Append(".\n");
        }
        else
        {
            if (sub is not null)
                sb.Append("    substitution_score = ").Append(Real(sub.Value)).Append(",\n");
            if (gap is not null)
                sb.Append("    gap_score = ").Append(Real(gap.Value)).Append(".\n");
        }

        string L(long v) => v.ToString(CultureInfo.InvariantCulture);
        sb.Append("    aligned = ").Append(L(Aligned)).Append(":\n");
        sb.Append("        identities = ").Append(L(Identities)).Append(",\n");
        if (Positives is not null)
            sb.Append("        positives = ").Append(L(Positives.Value)).Append(",\n");
        sb.Append("        mismatches = ").Append(L(Mismatches)).Append(".\n");
        sb.Append("    gaps = ").Append(L(Gaps)).Append(":\n");
        sb.Append("        left_gaps = ").Append(L(LeftGaps)).Append(":\n");
        sb.Append("            left_insertions = ").Append(L(LeftInsertions)).Append(":\n");
        sb.Append("                open_left_insertions = ").Append(L(OpenLeftInsertions)).Append(",\n");
        sb.Append("                extend_left_insertions = ").Append(L(ExtendLeftInsertions)).Append(";\n");
        sb.Append("            left_deletions = ").Append(L(LeftDeletions)).Append(":\n");
        sb.Append("                open_left_deletions = ").Append(L(OpenLeftDeletions)).Append(",\n");
        sb.Append("                extend_left_deletions = ").Append(L(ExtendLeftDeletions)).Append(";\n");
        sb.Append("        internal_gaps = ").Append(L(InternalGaps)).Append(":\n");
        sb.Append("            internal_insertions = ").Append(L(InternalInsertions)).Append(":\n");
        sb.Append("                open_internal_insertions = ").Append(L(OpenInternalInsertions)).Append(",\n");
        sb.Append("                extend_internal_insertions = ").Append(L(ExtendInternalInsertions)).Append(";\n");
        sb.Append("            internal_deletions = ").Append(L(InternalDeletions)).Append(":\n");
        sb.Append("                open_internal_deletions = ").Append(L(OpenInternalDeletions)).Append(",\n");
        sb.Append("                extend_internal_deletions = ").Append(L(ExtendInternalDeletions)).Append(";\n");
        sb.Append("        right_gaps = ").Append(L(RightGaps)).Append(":\n");
        sb.Append("            right_insertions = ").Append(L(RightInsertions)).Append(":\n");
        sb.Append("                open_right_insertions = ").Append(L(OpenRightInsertions)).Append(",\n");
        sb.Append("                extend_right_insertions = ").Append(L(ExtendRightInsertions)).Append(";\n");
        sb.Append("            right_deletions = ").Append(L(RightDeletions)).Append(":\n");
        sb.Append("                open_right_deletions = ").Append(L(OpenRightDeletions)).Append(",\n");
        sb.Append("                extend_right_deletions = ").Append(L(ExtendRightDeletions)).Append(".\n");
        return sb.ToString();
    }
}

public static partial class AlignmentOperations
{
    private const int Horizontal = 1, Vertical = 2, Diagonal = 4;

    /// <summary><c>alignment.counts()</c>: the counts without a scoring (identities and mismatches by letter equality; no scores).</summary>
    public static AlignmentCounts Counts(BioCoordAlignment alignment) => Counts(alignment, null, null, null);

    /// <summary><c>alignment.counts(wildcard)</c>: letters equal to <paramref name="wildcard"/> are neither identities nor mismatches.</summary>
    public static AlignmentCounts Counts(BioCoordAlignment alignment, char wildcard) => Counts(alignment, null, null, wildcard);

    /// <summary><c>alignment.counts(substitution_matrix)</c>: identities, mismatches and positives with the scores of the matrix (letters must be in its alphabet).</summary>
    public static AlignmentCounts Counts(BioCoordAlignment alignment, SubstitutionMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        return Counts(alignment, null, matrix, null);
    }

    /// <summary><c>alignment.counts(aligner)</c>: the counts together with the gap, substitution and total score under the scoring of the aligner (match / mismatch or matrix, wildcard, gap scores, gap-score functions).</summary>
    public static AlignmentCounts Counts(BioCoordAlignment alignment, PairwiseAligner aligner)
    {
        ArgumentNullException.ThrowIfNull(aligner);
        return Counts(alignment, aligner, aligner.SubstitutionMatrix, aligner.Wildcard);
    }

    private static AlignmentCounts Counts(BioCoordAlignment alignment, PairwiseAligner? aligner, SubstitutionMatrix? matrix, char? wildcard)
    {
        RequireCoordinates(alignment);
        int n = alignment.Sequences.Count;
        var coordinates = alignment.Coordinates.Select(r => (long[])r.Clone()).ToArray();
        var sequences = alignment.Sequences.ToArray();
        var strands = new bool[n];
        var steps = StepsOf(coordinates);
        var alignedFlags = AlignedBlocks(steps);
        for (int i = 0; i < n; i++)
        {
            int pos = 0, neg = 0;
            for (int k = 0; k < steps[i].Length; k++)
            {
                if (!alignedFlags[k])
                    continue;
                if (steps[i][k] > 0)
                    pos++;
                else if (steps[i][k] < 0)
                    neg++;
            }

            if (pos < neg)
            {
                long len = sequences[i].Length;
                sequences[i] = AlignmentCoordinateFormats.ReverseComplement(sequences[i]);
                for (int k = 0; k < coordinates[i].Length; k++)
                    coordinates[i][k] = len - coordinates[i][k];
                strands[i] = true;
            }
        }

        string?[] texts = sequences.Select(s => s.Parts is null ? null : Contents(s)).ToArray();
        long[] gapCount = new long[12];
        long aligned = 0, identities = 0, mismatches = 0, positives = matrix is null ? -1 : 0;
        double gapScore = 0, substitutionScore = 0;
        int path = 0;
        int columns = coordinates.Length == 0 ? 0 : coordinates[0].Length;
        for (int jA = 0; jA < n; jA++)
        {
            for (int jB = jA + 1; jB < n; jB++)
            {
                long leftA = coordinates[jA][0], leftB = coordinates[jB][0];
                long rightA = coordinates[jA][columns - 1], rightB = coordinates[jB][columns - 1];
                long startA = leftA, startB = leftB;
                for (int k = 1; k < columns; k++)
                {
                    long endA = coordinates[jA][k], endB = coordinates[jB][k];
                    if (startA == endA && startB == endB)
                    {
                        // nothing
                    }
                    else if (startA == endA)
                    {
                        AddGaps(ref path, Horizontal, strands[jA], aligner?.InsertionScoreFunction, startA, endA, leftA, rightA, endB - startB, gapCount, 0, ref gapScore);
                    }
                    else if (startB == endB)
                    {
                        AddGaps(ref path, Vertical, strands[jB], aligner?.DeletionScoreFunction, startB, endB, leftB, rightB, endA - startA, gapCount, 1, ref gapScore);
                    }
                    else if (texts[jA] is null || texts[jB] is null)
                    {
                        path = Diagonal;
                        aligned += endA - startA;
                    }
                    else
                    {
                        path = Diagonal;
                        aligned += endA - startA;
                        string a = texts[jA]!, b = texts[jB]!;
                        for (long lA = startA, lB = startB; lA < endA && lB < endB; lA++, lB++)
                        {
                            char cA = a[(int)lA], cB = b[(int)lB];
                            if (matrix is null)
                            {
                                if (wildcard is { } w && (cA == w || cB == w))
                                    continue;
                                if (cA == cB)
                                    identities++;
                                else
                                    mismatches++;
                            }
                            else
                            {
                                if (!matrix.TryGetValue(cA, cB, out decimal value))
                                    throw Bad("sequence contains letters not in the alphabet");
                                if (cA == cB)
                                    identities++;
                                else
                                    mismatches++;
                                if (value > 0)
                                    positives++;
                                substitutionScore += (double)value;
                            }
                        }
                    }

                    startA = endA;
                    startB = endB;
                }
            }
        }

        double? substitution = null, gap = null;
        if (matrix is not null)
        {
            substitution = substitutionScore;
        }
        else if (aligner is not null)
        {
            substitution = aligner.MatchScore * identities + aligner.MismatchScore * mismatches;
        }

        if (aligner is not null)
        {
            double g = gapScore;
            if (aligner.InsertionScoreFunction is null)
                g += gapCount[0] * aligner.OpenLeftInsertionScore + gapCount[1] * aligner.ExtendLeftInsertionScore
                   + gapCount[2] * aligner.OpenInternalInsertionScore + gapCount[3] * aligner.ExtendInternalInsertionScore
                   + gapCount[4] * aligner.OpenRightInsertionScore + gapCount[5] * aligner.ExtendRightInsertionScore;
            if (aligner.DeletionScoreFunction is null)
                g += gapCount[6] * aligner.OpenLeftDeletionScore + gapCount[7] * aligner.ExtendLeftDeletionScore
                   + gapCount[8] * aligner.OpenInternalDeletionScore + gapCount[9] * aligner.ExtendInternalDeletionScore
                   + gapCount[10] * aligner.OpenRightDeletionScore + gapCount[11] * aligner.ExtendRightDeletionScore;
            gap = g;
        }

        return new AlignmentCounts
        {
            OpenLeftInsertions = gapCount[0],
            ExtendLeftInsertions = gapCount[1],
            OpenInternalInsertions = gapCount[2],
            ExtendInternalInsertions = gapCount[3],
            OpenRightInsertions = gapCount[4],
            ExtendRightInsertions = gapCount[5],
            OpenLeftDeletions = gapCount[6],
            ExtendLeftDeletions = gapCount[7],
            OpenInternalDeletions = gapCount[8],
            ExtendInternalDeletions = gapCount[9],
            OpenRightDeletions = gapCount[10],
            ExtendRightDeletions = gapCount[11],
            Aligned = aligned,
            Identities = identities,
            Mismatches = mismatches,
            Positives = matrix is null ? null : positives,
            SubstitutionScore = substitution,
            GapScore = gap,
        };
    }

    // add_gaps of _alignmentcounts.c; the counters are [open left, extend left, open internal, extend internal, open right, extend right] for insertions (offset 0) and deletions (offset 6)
    private static void AddGaps(ref int path, int direction, bool strand, Func<int, int, double>? scoreFunction, long start, long end, long left, long right, long gapSize, long[] counters, int kind, ref double gapScore)
    {
        int o = kind * 6;
        if (path == direction)
        {
            if (start == left)
                counters[o + 1] += gapSize;
            else if (end == right)
                counters[o + 5] += gapSize;
            else
                counters[o + 3] += gapSize;
            return;
        }

        if (start == left)
        {
            counters[o]++;
            counters[o + 1] += gapSize - 1;
        }
        else if (end == right)
        {
            counters[o + 4]++;
            counters[o + 5] += gapSize - 1;
        }
        else
        {
            counters[o + 2]++;
            counters[o + 3] += gapSize - 1;
        }

        if (scoreFunction is not null)
            gapScore += strand ? scoreFunction((int)(right - start), (int)gapSize) : scoreFunction((int)start, (int)gapSize);
        path = direction;
    }
}
