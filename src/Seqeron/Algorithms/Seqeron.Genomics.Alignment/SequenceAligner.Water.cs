using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

// EMBOSS 6.6 water (audit A2-28, B09 F122): nucleus/embaln.c embAlignPathCalcSW + embAlignWalkSWMatrix with the program's arithmetic (float path matrix, double running maxima) and its walk-back
// (the cell of the best score is the LAST one scanned within U_FEPS of the maximum; a gap is followed back to the cell that explains the score within 0.01), so the rows are the rows of the program.
public static partial class SequenceAligner
{
    /// <summary>
    /// EMBOSS <c>water</c>: local alignment (Smith-Waterman, Gotoh's affine gaps) as the program computes and reports it: a gap of length k costs <c>gapOpen + (k - 1) * gapExtend</c>
    /// (positive penalties, <c>-gapopen</c> / <c>-gapextend</c>, rounded to 16 mantissa bits as EMBOSS does), scores are single-precision. The rows (ties included), the score and the
    /// coordinates are those of the installed water 6.6.0 (the program needs sequences of at least two residues: it reports its path-matrix maximum, which is undefined for a one-row matrix;
    /// shorter input returns the empty alignment).
    /// </summary>
    /// <param name="sequence1">First sequence (<c>-asequence</c>).</param>
    /// <param name="sequence2">Second sequence (<c>-bsequence</c>).</param>
    /// <param name="matrix">Residue-pair scores (<c>-datafile</c>); default EDNAFULL for nucleotide sequences, BLOSUM62 (EBLOSUM62) otherwise.</param>
    /// <param name="gapOpen">Cost of a gap of length one (<c>-gapopen</c>, default 10).</param>
    /// <param name="gapExtend">Cost of each further gap position (<c>-gapextend</c>, default 0.5).</param>
    /// <returns>The local alignment: its rows, <see cref="AlignmentResult.StartPosition1"/> / <see cref="AlignmentResult.StartPosition2"/> (0-based first aligned residues), the ends, and the single-precision score.</returns>
    public static AlignmentResult WaterAlign(string sequence1, string sequence2, SubstitutionMatrix? matrix = null, double gapOpen = 10, double gapExtend = 0.5)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        RequireEmbossRange(gapOpen, 1, 100, nameof(gapOpen), "-gapopen", "water");
        RequireEmbossRange(gapExtend, 0, 10, nameof(gapExtend), "-gapextend", "water");
        if (sequence1.Length < 2 || sequence2.Length < 2)
            return new AlignmentResult("", "", 0, AlignmentType.Local, 0, 0, 0, 0);
        matrix ??= IsNucleotideText(sequence1.ToUpperInvariant()) && IsNucleotideText(sequence2.ToUpperInvariant()) ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62;
        var (score, r1, r2, start1, start2) = WaterCore(sequence1, sequence2, matrix, RoundFloat((float)gapOpen, 8), RoundFloat((float)gapExtend, 8));
        int len1 = r1.Count(c => c != '-');
        int len2 = r2.Count(c => c != '-');
        return new AlignmentResult(r1, r2, (int)Math.Round(score, MidpointRounding.AwayFromZero), AlignmentType.Local, start1, start2, start1 + len1, start2 + len2) { ExactScore = score };
    }

    private static (float Score, string Row1, string Row2, int Start1, int Start2) WaterCore(string s1, string s2, SubstitutionMatrix matrix, float gapopen, float gapextend)
    {
        int lena = s1.Length;
        int lenb = s2.Length;
        RequireCells((long)lena * lenb, WaterMaxCells, nameof(WaterMaxCells), "EMBOSS water", "keeps a float path matrix and an int compass matrix (8 bytes per cell)", $"{lena:N0} x {lenb:N0} residues = {(long)lena * lenb:N0} cells",
            "The program has no fixed limit: it dies with \"Sequences too big. Try 'matcher' or 'supermatcher'\" when the allocation fails", "use LocalAlignLinearSpace / Supermatcher for sequences of this length");
        var table = new EmbossScoreTable(matrix);
        int[] ca = table.Codes(s1);
        int[] cb = table.Codes(s2);
        float Sub(int x, int y) => table.Score(x, y);

        long n = (long)lena * lenb;
        var path = new float[n];
        var compass = new int[n];
        var maxa = new double[lena];
        float ret = -float.MaxValue;

        for (int i = 0; i < lena; i++)
        {
            double result = Sub(ca[i], cb[0]);
            double fnew;
            if (i == 0)
            {
                fnew = 0.0;
            }
            else
            {
                float t = path[(long)(i - 1) * lenb] - (compass[(long)(i - 1) * lenb] == 2 ? gapextend : gapopen);
                fnew = t;
            }

            if (result > fnew && result > 0)
            {
                path[(long)i * lenb] = (float)result;
                compass[(long)i * lenb] = 0;
            }
            else if (fnew > 0)
            {
                path[(long)i * lenb] = (float)fnew;
                compass[(long)i * lenb] = 2;
            }
            else
            {
                path[(long)i * lenb] = 0f;
                compass[(long)i * lenb] = 0;
            }

            float mx = path[(long)i * lenb] - gapopen;
            if (i > 0 && compass[(long)(i - 1) * lenb] == 2)
                mx = path[(long)i * lenb] - gapextend;
            maxa[i] = mx;
        }

        for (int j = 0; j < lenb; j++)
        {
            double result = Sub(ca[0], cb[j]);
            double fnew;
            if (j == 0)
            {
                fnew = 0.0;
            }
            else
            {
                float t = path[j - 1] - (compass[j - 1] == 1 ? gapextend : gapopen);
                fnew = t;
            }

            if (result > fnew && result > 0)
            {
                path[j] = (float)result;
                compass[j] = 0;
            }
            else if (fnew > 0)
            {
                path[j] = (float)fnew;
                compass[j] = 1;
            }
            else
            {
                path[j] = 0f;
                compass[j] = 0;
            }
        }

        for (int xpos = 1; xpos != lenb; xpos++)
        {
            AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per column of the matrix
            float bx0 = path[xpos] - gapopen - gapextend;
            double bx = bx0;
            for (int ypos = 1; ypos < lena; ypos++)
            {
                double match = Sub(ca[ypos], cb[xpos]);
                double mscore = path[(long)(ypos - 1) * lenb + xpos - 1] + match;
                long cell = ((long)ypos * lenb) + xpos;
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
                    compass[cell] = 1;
                }

                bx -= gapextend;
                fnew = path[cell - lenb];
                fnew -= gapopen;
                if (fnew > bx)
                    bx = fnew;
                if (bx > mscore)
                {
                    mscore = bx;
                    path[cell] = (float)mscore;
                    compass[cell] = 2;
                }

                if (mscore > ret)
                    ret = (float)mscore;
                double rr = path[cell];
                if (rr < 0.0)
                    path[cell] = 0f;
            }
        }

        // embAlignWalkSWMatrix
        double pmax = -float.MaxValue;
        long k = n - 1;
        long xp = 0, yp = 0;
        for (int i = lena - 1; i >= 0; --i)
        {
            for (int j = lenb - 1; j >= 0; --j)
            {
                float here = path[k];
                k--;
                bool eq = (pmax - NeedleFloatEpsilon) < here && here < (pmax + NeedleFloatEpsilon);
                if (here > pmax || eq)
                {
                    pmax = here;
                    xp = j;
                    yp = i;
                }
            }
        }

        var mrow = new StringBuilder();
        var nrow = new StringBuilder();
        const double errbounds = 0.01;
        while (xp >= 0 && yp >= 0)
        {
            long cur = (yp * lenb) + xp;
            if (compass[cur] == 0)
            {
                mrow.Append(s1[(int)yp--]);
                nrow.Append(s2[(int)xp--]);
                if (yp >= 0 && xp >= 0 && path[(yp * lenb) + xp] <= 0f)
                    break;
            }
            else if (compass[cur] == 1)
            {
                double score = path[cur];
                long gapcnt = 0;
                long ix = xp - 1;
                float bimble;
                while (true)
                {
                    bimble = path[(yp * lenb) + ix] - gapopen - (gapcnt * gapextend);
                    if (ix == 0 || Math.Abs(score - (double)bimble) < errbounds)
                        break;
                    --ix;
                    ++gapcnt;
                }

                if (bimble <= 0.0)
                    break;
                for (long ic = 0; ic <= gapcnt; ++ic)
                {
                    mrow.Append('-');
                    nrow.Append(s2[(int)xp--]);
                }
            }
            else
            {
                double score = path[cur];
                long gapcnt = 0;
                long iy = yp - 1;
                float bimble;
                while (true)
                {
                    bimble = path[(iy * lenb) + xp] - gapopen - (gapcnt * gapextend);
                    if (iy == 0 || Math.Abs(score - (double)bimble) < errbounds)
                        break;
                    --iy;
                    if (iy < 0)
                        throw new InvalidOperationException("SW: Error walking down");
                    ++gapcnt;
                }

                if (bimble <= 0.0)
                    break;
                for (long ic = 0; ic <= gapcnt; ++ic)
                {
                    mrow.Append(s1[(int)yp--]);
                    nrow.Append('-');
                }
            }
        }

        return (ret, Reverse(mrow), Reverse(nrow), (int)(yp + 1), (int)(xp + 1));
    }
}
