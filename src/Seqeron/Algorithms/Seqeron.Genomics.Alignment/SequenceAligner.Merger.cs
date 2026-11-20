using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>A disagreement of <see cref="SequenceAligner.Merger"/>: a column where the two sequences differ (or one has a gap) and the base that went into the merged sequence.</summary>
/// <param name="PositionA">1-based position in the first sequence (the next residue of sequence A at the column).</param>
/// <param name="BaseA">The base of sequence A at the column ('.' for a gap, as the program prints it).</param>
/// <param name="PositionB">1-based position in the second sequence.</param>
/// <param name="BaseB">The base of sequence B ('.' for a gap).</param>
/// <param name="Using">The base taken for the merged sequence ('.' when a gap was chosen).</param>
public readonly record struct MergerConflict(int PositionA, char BaseA, int PositionB, char BaseB, char Using);

/// <summary>Result of <see cref="SequenceAligner.Merger"/> (EMBOSS merger).</summary>
/// <param name="MergedSequence">The merged sequence; letters taken from a conflict are upper case, everything else is lower case, as the program writes it.</param>
/// <param name="Score">The alignment score (single precision, as <c>embAlignPathCalc</c> returns it).</param>
/// <param name="StartA">0-based number of residues of sequence A before the aligned region (the program's <c>start1</c>).</param>
/// <param name="StartB">The same for sequence B (<c>start2</c>).</param>
/// <param name="AlignedA">The aligned region of sequence A with '.' for gaps (the program's walk string).</param>
/// <param name="AlignedB">The aligned region of sequence B.</param>
/// <param name="Conflicts">The conflict lines of the report, in order.</param>
/// <param name="Warnings">The warning lines of the report ("The region of alignment ...", "There is an equal amount of unaligned sequence ...").</param>
public sealed record MergerResult(
    string MergedSequence, float Score, int StartA, int StartB, string AlignedA, string AlignedB,
    IReadOnlyList<MergerConflict> Conflicts, IReadOnlyList<string> Warnings);

public static partial class SequenceAligner
{
    #region merger (merge two overlapping sequences)

    /// <summary>
    /// EMBOSS <c>merger</c> (6.6.0, <c>emboss/merger.c</c>): merges two overlapping nucleotide sequences into one. A global alignment with free end gaps
    /// (<c>embAlignPathCalc</c> + <c>embAlignWalkNWMatrix</c>: Gotoh-style row / column maxima in single-precision path matrices, a gap of length k costing
    /// <c>gapOpen + (k - 1) * gapExtend</c>, the alignment starting and ending on the matrix border) is walked column by column: the unaligned head of the longer-leading
    /// sequence, the aligned columns and the unaligned tail of the longer-trailing sequence make the merged sequence.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Column rule: identical bases are copied; an <c>N</c> in one sequence yields the other sequence's base; a mismatch or a gap takes the base of the sequence with the better
    /// <b>quality</b> around the position (<c>merger_quality</c>: in a window of w = 5, then 10, then 20 positions on both sides, 2 points for A/C/G/T/U, 1 for another IUPAC code,
    /// 0 for N and for positions off the end of the sequence, divided by 2w + 1; the first window size at which the two differ decides, equal quality picks sequence A) and records
    /// the conflict. If both sequences have unaligned residues at the start (end) of the alignment only one of them is used (the longer one; the first when equal) and a warning
    /// is recorded. Both sequences are lower-cased first; the chosen base of a conflict (and the base kept next to an N) is upper-cased in the result.
    /// </para>
    /// <para>
    /// The scoring matrix is EDNAFULL for nucleotide sequences (the program's default <c>-datafile</c>) or any <see cref="SubstitutionMatrix"/>. Verified against merger 6.6.0 on 400 random
    /// overlapping read pairs (merged sequence, conflict list, warnings; docs/Evidence/ALIGN-SEMI-001-Evidence.md).
    /// </para>
    /// </remarks>
    /// <param name="sequenceA">First sequence (<c>-asequence</c>).</param>
    /// <param name="sequenceB">Second sequence (<c>-bsequence</c>).</param>
    /// <param name="matrix">Scoring matrix (default EDNAFULL).</param>
    /// <param name="gapOpen">Gap opening penalty (<c>-gapopen</c>, 0..100, default 50).</param>
    /// <param name="gapExtend">Gap extension penalty (<c>-gapextend</c>, 0..10, default 5).</param>
    /// <param name="cancellationToken">Cancellation token, polled per column of the path matrix and per step of the walk (A4-1).</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentException">A sequence is empty.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A penalty is outside the program's range.</exception>
    public static MergerResult Merger(string sequenceA, string sequenceB, SubstitutionMatrix? matrix = null, double gapOpen = 50, double gapExtend = 5, CancellationToken cancellationToken = default)
    {
        using var scope = AlignmentCancellation.Enter(cancellationToken);
        ArgumentNullException.ThrowIfNull(sequenceA);
        ArgumentNullException.ThrowIfNull(sequenceB);
        sequenceA = EmbossSequenceReader.Clean(sequenceA);
        sequenceB = EmbossSequenceReader.Clean(sequenceB);
        RequireEmbossRange(gapOpen, 0, 100, nameof(gapOpen), "-gapopen", "merger");
        RequireEmbossRange(gapExtend, 0, 10, nameof(gapExtend), "-gapextend", "merger");
        if (sequenceA.Length == 0 || sequenceB.Length == 0)
            throw new ArgumentException("EMBOSS merger cannot merge an empty sequence.", sequenceA.Length == 0 ? nameof(sequenceA) : nameof(sequenceB));
        string a = sequenceA.ToLowerInvariant();
        string b = sequenceB.ToLowerInvariant();
        matrix ??= SubstitutionMatrix.EdnaFull;
        RequireCells((long)a.Length * b.Length, NeedleMaxCells, nameof(NeedleMaxCells), "EMBOSS merger", "keeps a float path matrix and an int compass matrix (8 bytes per cell)",
            $"{a.Length:N0} x {b.Length:N0} residues", "The program dies with \"Sequences too big\" when the allocation fails", "use megamerger-style word matching for sequences of this length");

        float go = RoundFloat((float)gapOpen, 8), ge = RoundFloat((float)gapExtend, 8);
        var (path, compass, score) = MergerPathCalc(a, b, go, ge, matrix);
        var (m, n, start1, start2) = MergerWalk(path, compass, a, b, go, ge);
        return MergerMerge(a, b, m, n, start1, start2, score);
    }

    // embAlignPathCalc: path(i, j) = best score of an alignment ending at a[i] / b[j]; maxa[i] / maxb = best open gap along the row / column (double), path (float).
    private static (float[] Path, int[] Compass, float Score) MergerPathCalc(string a, string b, float gapopen, float gapextend, SubstitutionMatrix matrix)
    {
        const int Left = 1, Down = 2;
        int lena = a.Length, lenb = b.Length;
        var table = new EmbossScoreTable(matrix);
        int[] ca = table.Codes(a), cb = table.Codes(b);
        var path = new float[(long)lena * lenb];
        var compass = new int[(long)lena * lenb];
        var maxa = new double[lena];

        for (int i = 0; i < lena; i++)
        {
            path[(long)i * lenb] = table.Score(ca[i], cb[0]);
            maxa[i] = path[(long)i * lenb] - (double)gapopen;
        }

        for (int j = 0; j < lenb; j++)
            path[j] = table.Score(ca[0], cb[j]);

        for (int xpos = 1; xpos < lenb; xpos++)
        {
            AlignmentCancellation.ThrowIfRequested();
            double maxb = path[xpos] - (double)gapopen;
            for (int ypos = 1; ypos < lena; ypos++)
            {
                long cell = (long)ypos * lenb + xpos;
                double match = table.Score(ca[ypos], cb[xpos]);
                double mscore = path[cell - lenb - 1] + match;
                compass[cell] = 0;
                path[cell] = (float)mscore;

                maxa[ypos] -= gapextend;
                double fnew = path[cell - 1];
                fnew -= gapopen;
                if (fnew > maxa[ypos])
                    maxa[ypos] = fnew;
                if (maxa[ypos] > mscore)
                {
                    mscore = maxa[ypos];
                    path[cell] = (float)mscore;
                    compass[cell] = Left;
                }

                maxb -= gapextend;
                fnew = path[cell - lenb];
                fnew -= gapopen;
                if (fnew > maxb)
                    maxb = fnew;
                if (maxb > mscore)
                {
                    mscore = maxb;
                    path[cell] = (float)mscore;
                    compass[cell] = Down;
                }
            }
        }

        float ret = -float.MaxValue;
        for (int i = 0; i < lenb; i++)
            if (path[(long)(lena - 1) * lenb + i] > ret) ret = path[(long)(lena - 1) * lenb + i];
        for (int j = 0; j < lena; j++)
            if (path[(long)j * lenb + lenb - 1] > ret) ret = path[(long)j * lenb + lenb - 1];
        return (path, compass, ret);
    }

    // embAlignWalkNWMatrix: from the best border cell (last row, then last column; ties: the later cell of the row scan, the first of the column scan) back to the border.
    private static (string M, string N, int Start1, int Start2) MergerWalk(float[] path, int[] compass, string p, string q, float gapopen, float gapextend)
    {
        const int Left = 1, Down = 2;
        int lena = p.Length, lenb = q.Length;
        const double errbounds = 0.01;
        double pmax = -1.0 * int.MaxValue;
        int xpos = 0, ypos = 0;
        for (int i = 0; i < lenb; i++)
            if (path[(long)(lena - 1) * lenb + i] >= pmax)
            {
                pmax = path[(long)(lena - 1) * lenb + i];
                xpos = i;
                ypos = lena - 1;
            }

        for (int j = 0; j < lena; j++)
            if (path[(long)j * lenb + lenb - 1] > pmax)
            {
                pmax = path[(long)j * lenb + lenb - 1];
                xpos = lenb - 1;
                ypos = j;
            }

        var m = new StringBuilder();
        var n = new StringBuilder();
        int steps = 0;
        while (xpos >= 0 && ypos >= 0)
        {
            if ((++steps & 1023) == 0)
                AlignmentCancellation.ThrowIfRequested();
            int dir = compass[(long)ypos * lenb + xpos];
            if (dir == 0)
            {
                m.Append(p[ypos--]);
                n.Append(q[xpos--]);
            }
            else if (dir == Left)
            {
                double score = path[(long)ypos * lenb + xpos];
                int gapcnt = 0, ix = xpos - 1;
                while (true)
                {
                    double bimble = path[(long)ypos * lenb + ix] - (double)gapopen - gapcnt * (double)gapextend;
                    if (ix == 0 || Math.Abs(score - bimble) < errbounds)
                        break;
                    ix--;
                    gapcnt++;
                }

                for (int ic = 0; ic <= gapcnt; ic++)
                {
                    m.Append('.');
                    n.Append(q[xpos--]);
                }
            }
            else if (dir == Down)
            {
                double score = path[(long)ypos * lenb + xpos];
                int gapcnt = 0, iy = ypos - 1;
                while (true)
                {
                    double bimble = path[(long)iy * lenb + xpos] - (double)gapopen - gapcnt * (double)gapextend;
                    if (iy == 0 || Math.Abs(score - bimble) < errbounds)
                        break;
                    iy--;
                    gapcnt++;
                }

                for (int ic = 0; ic <= gapcnt; ic++)
                {
                    m.Append(p[ypos--]);
                    n.Append('.');
                }
            }
            else
            {
                throw new InvalidOperationException("Walk Error in NW");
            }
        }

        var ma = m.ToString().ToCharArray();
        var na = n.ToString().ToCharArray();
        Array.Reverse(ma);
        Array.Reverse(na);
        return (new string(ma), new string(na), ypos + 1, xpos + 1);
    }

    // merger_Merge + merger_bestquality + merger_quality.
    private static MergerResult MergerMerge(string a, string b, string m, string n, int start1, int start2, float score)
    {
        var ms = new StringBuilder();
        var conflicts = new List<MergerConflict>();
        var warnings = new List<string>();
        var p = m.ToCharArray();
        var q = n.ToCharArray();
        int olen = p.Length;

        if (start1 > start2)
        {
            ms.Append(a, 0, start1);
            if (start2 != 0)
            {
                warnings.Add($"The region of alignment only starts at position {start2 + 1} of sequence B");
                warnings.Add("Only the sequence of A is being used before this point");
            }
        }
        else if (start2 > start1)
        {
            ms.Append(b, 0, start2);
            if (start1 != 0)
            {
                warnings.Add($"The region of alignment only starts at position {start1 + 1} of sequence A");
                warnings.Add("Only the sequence of B is being used before this point");
            }
        }
        else if (start1 != 0 && start2 != 0)
        {
            warnings.Add($"There is an equal amount of unaligned sequence ({start1}) at the start of the sequences.");
            warnings.Add("Sequence A is being arbitrarily chosen for the merged sequence");
            ms.Append(a, 0, start1);
        }

        static bool IsGap(char c) => c is '.' or ' ' or '-';
        int apos = start1, bpos = start2;
        for (int i = 0; i < olen; i++)
        {
            if (IsGap(p[i]) || IsGap(q[i]))
            {
                if (MergerBestQuality(a, b, apos, bpos))
                {
                    p[i] = char.ToUpperInvariant(p[i]);
                    if (!IsGap(p[i])) ms.Append(p[i]);
                    conflicts.Add(new MergerConflict(apos + 1, p[i], bpos + 1, q[i], p[i]));
                }
                else
                {
                    q[i] = char.ToUpperInvariant(q[i]);
                    if (!IsGap(q[i])) ms.Append(q[i]);
                    conflicts.Add(new MergerConflict(apos + 1, p[i], bpos + 1, q[i], q[i]));
                }
            }
            else if (p[i] is 'n' or 'N')
            {
                q[i] = char.ToUpperInvariant(q[i]);
                ms.Append(q[i]);
            }
            else if (q[i] is 'n' or 'N')
            {
                p[i] = char.ToUpperInvariant(p[i]);
                ms.Append(p[i]);
            }
            else if (p[i] != q[i])
            {
                if (MergerBestQuality(a, b, apos, bpos))
                {
                    p[i] = char.ToUpperInvariant(p[i]);
                    ms.Append(p[i]);
                    conflicts.Add(new MergerConflict(apos + 1, p[i], bpos + 1, q[i], p[i]));
                }
                else
                {
                    q[i] = char.ToUpperInvariant(q[i]);
                    ms.Append(q[i]);
                    conflicts.Add(new MergerConflict(apos + 1, p[i], bpos + 1, q[i], q[i]));
                }
            }
            else
            {
                ms.Append(p[i]);
            }

            if (!IsGap(p[i])) apos++;
            if (!IsGap(q[i])) bpos++;
        }

        int alen = a.Length - apos, blen = b.Length - bpos;
        if (alen > blen)
        {
            ms.Append(a, apos, alen);
            if (blen != 0)
            {
                warnings.Add($"The region of alignment ends at position {bpos + 1} of sequence B");
                warnings.Add("Only the sequence of A is being used after this point");
            }
        }

        if (blen > alen)
        {
            ms.Append(b, bpos, blen);
            if (alen != 0)
            {
                warnings.Add($"The region of alignment ends at position {apos + 1} of sequence A");
                warnings.Add("Only the sequence of B is being used after this point");
            }
        }
        else if (alen != 0 && blen != 0)
        {
            // merger.c: `if(blen > alen) {...} else if(alen && blen) {...}` - this branch also runs after the alen > blen block, so such a tail is appended twice (a program quirk, reproduced)
            warnings.Add($"There is an equal amount of unaligned sequence ({alen}) at the end of the sequences.");
            warnings.Add("Sequence A is being arbitrarily chosen for the merged sequence");
            ms.Append(a, apos, alen);
        }

        return new MergerResult(ms.ToString(), score, start1, start2, m, n, conflicts, warnings);
    }

    private static bool MergerBestQuality(string a, string b, int apos, int bpos)
    {
        float qa = MergerQuality(a, apos, 5), qb = MergerQuality(b, bpos, 5);
        if (qa == qb)
        {
            qa = MergerQuality(a, apos, 10);
            qb = MergerQuality(b, bpos, 10);
        }

        if (qa == qb)
        {
            qa = MergerQuality(a, apos, 20);
            qb = MergerQuality(b, bpos, 20);
        }

        return qa >= qb;
    }

    private static float MergerQuality(string seq, int pos, int window)
    {
        int value = 0;
        static int Points(char c)
        {
            if ("aAcCgGtTuU".Contains(c, StringComparison.Ordinal)) return 2;
            return "mrwsykvhdbMRWSYKVHDB".Contains(c, StringComparison.Ordinal) ? 1 : 0;
        }

        for (int i = pos; i < pos + window && i < seq.Length; i++)
            value += Points(seq[i]);
        for (int j = pos - 1; j > pos - window && j >= 0; j--)
            value += Points(seq[j]);
        return (float)((double)value / ((double)window * 2.0 + 1.0));
    }

    #endregion
}
