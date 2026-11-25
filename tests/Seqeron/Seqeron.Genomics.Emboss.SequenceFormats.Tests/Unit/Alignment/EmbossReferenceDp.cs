// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// Independent, score-only dynamic programs of the two EMBOSS 6.6.0 recurrences (audit A3-14b, B09 F143), written from the C text of <c>nucleus/embaln.c</c> in plain doubles
/// (no traceback, no tie rules, no float rounding: the penalties they are used with are multiples of 1/8 and the matrix scores are integers, so every value is exact): they say
/// what the PROGRAM's recurrence computes, which differs from the textbook Gotoh optimum in exactly the places the C code differs from it, so that a difference between the port
/// and the plain optimum is explained by a recomputation, not by a tolerance.
/// <list type="bullet">
/// <item><c>embAlignPathCalcWithEndGapPenalties</c> (needle): the last row / column use only the match state as the source of an end gap (a gap of the other direction cannot be followed by an end gap), the
/// first row / column are initialised with the end-gap penalties.</item>
/// <item><c>embAlignPathCalcSW</c> (water): the first row and column are initialised separately and the running maximum <c>ret</c> is taken only over the cells with both indices at least 1, so an
/// alignment that ENDS in the first row or column is never reported; the maximum is taken before the floor at zero (all-negative matrices give a negative score).</item>
/// </list>
/// </summary>
internal static class EmbossReferenceDp
{
    public static double Needle(string s1, string s2, SubstitutionMatrix matrix, double gapOpen, double gapExtend, bool endWeight, double endOpen, double endExtend)
    {
        int lena = s1.Length, lenb = s2.Length;
        double Sub(int y, int x) => matrix.TryGetScore(char.ToUpperInvariant(s1[y]), char.ToUpperInvariant(s2[x]), out int sc) ? sc : 0;   // a residue outside the matrix scores 0 (EMBOSS 'bad character')
        if (!endWeight)
        {
            endOpen = 0;
            endExtend = 0;
        }

        var m = new double[lena, lenb];
        var ix = new double[lena, lenb];
        var iy = new double[lena, lenb];
        m[0, 0] = Sub(0, 0);
        ix[0, 0] = iy[0, 0] = -endOpen - gapOpen;
        for (int y = 1; y < lena; y++)
        {
            iy[y, 0] = Math.Max(m[y - 1, 0] - gapOpen, iy[y - 1, 0] - gapExtend);
            m[y, 0] = Sub(y, 0) - (endOpen + ((y - 1) * endExtend));
            ix[y, 0] = -endOpen - (y * endExtend) - gapOpen;
        }

        ix[lena - 1, 0] = ix[lena - 1, 0] - endOpen + gapOpen;
        for (int x = 1; x < lenb; x++)
        {
            ix[0, x] = Math.Max(m[0, x - 1] - gapOpen, ix[0, x - 1] - gapExtend);
            m[0, x] = Sub(0, x) - (endOpen + ((x - 1) * endExtend));
            iy[0, x] = -endOpen - (x * endExtend) - gapOpen;
        }

        iy[0, lenb - 1] = iy[0, lenb - 1] - endOpen + gapOpen;
        for (int x = 1; x < lenb; x++)
        {
            for (int y = 1; y < lena; y++)
            {
                m[y, x] = Math.Max(m[y - 1, x - 1], Math.Max(ix[y - 1, x - 1], iy[y - 1, x - 1])) + Sub(y, x);
                double og, eg;
                if (x == lenb - 1)
                {
                    og = m[y - 1, x] - endOpen;
                    eg = iy[y - 1, x] - endExtend;
                }
                else
                {
                    og = Math.Max(m[y - 1, x], ix[y - 1, x]) - gapOpen;
                    eg = iy[y - 1, x] - gapExtend;
                }

                iy[y, x] = Math.Max(og, eg);
                if (y == lena - 1)
                {
                    og = m[y, x - 1] - endOpen;
                    eg = ix[y, x - 1] - endExtend;
                }
                else
                {
                    og = Math.Max(m[y, x - 1], iy[y, x - 1]) - gapOpen;
                    eg = ix[y, x - 1] - gapExtend;
                }

                ix[y, x] = Math.Max(og, eg);
            }
        }

        return Math.Max(m[lena - 1, lenb - 1], Math.Max(ix[lena - 1, lenb - 1], iy[lena - 1, lenb - 1]));
    }

    /// <param name="includeBorderCells">false = the program (maximum over the cells with both indices &gt;= 1); true = the maximum also over the first row and column (the textbook local optimum, floored at zero).</param>
    public static double Water(string s1, string s2, SubstitutionMatrix matrix, double gapOpen, double gapExtend, bool includeBorderCells = false)
    {
        int lena = s1.Length, lenb = s2.Length;
        double Sub(int y, int x) => matrix.TryGetScore(char.ToUpperInvariant(s1[y]), char.ToUpperInvariant(s2[x]), out int sc) ? sc : 0;   // a residue outside the matrix scores 0 (EMBOSS 'bad character')
        const int Diag = 0, Left = 1, Down = 2;
        var path = new double[lena, lenb];
        var compass = new int[lena, lenb];
        var maxa = new double[lena];
        for (int i = 0; i < lena; i++)
        {
            double result = Sub(i, 0);
            double fnew = i == 0 ? 0 : path[i - 1, 0] - (compass[i - 1, 0] == Down ? gapExtend : gapOpen);
            if (result > fnew && result > 0) { path[i, 0] = result; compass[i, 0] = Diag; }
            else if (fnew > 0) { path[i, 0] = fnew; compass[i, 0] = Down; }
            else { path[i, 0] = 0; compass[i, 0] = Diag; }
            maxa[i] = i == 0 ? path[i, 0] - gapOpen : path[i, 0] - (compass[i - 1, 0] == Down ? gapExtend : gapOpen);
        }

        for (int j = 0; j < lenb; j++)
        {
            double result = Sub(0, j);
            double fnew = j == 0 ? 0 : path[0, j - 1] - (compass[0, j - 1] == Left ? gapExtend : gapOpen);
            if (result > fnew && result > 0) { path[0, j] = result; compass[0, j] = Diag; }
            else if (fnew > 0) { path[0, j] = fnew; compass[0, j] = Left; }
            else { path[0, j] = 0; compass[0, j] = Diag; }
        }

        double ret = double.MinValue;
        for (int x = 1; x < lenb; x++)
        {
            double bx = path[0, x] - gapOpen - gapExtend;
            for (int y = 1; y < lena; y++)
            {
                double mscore = path[y - 1, x - 1] + Sub(y, x);
                path[y, x] = mscore;
                maxa[y] -= gapExtend;
                double fnew = path[y, x - 1] - gapOpen;
                if (fnew > maxa[y]) maxa[y] = fnew;
                if (maxa[y] > mscore) { mscore = maxa[y]; path[y, x] = mscore; }
                bx -= gapExtend;
                fnew = path[y - 1, x] - gapOpen;
                if (fnew > bx) bx = fnew;
                if (bx > mscore) { mscore = bx; path[y, x] = mscore; }
                if (mscore > ret) ret = mscore;
                if (path[y, x] < 0) path[y, x] = 0;
            }
        }

        if (includeBorderCells)
        {
            for (int i = 0; i < lena; i++)
                ret = Math.Max(ret, path[i, 0]);
            for (int j = 0; j < lenb; j++)
                ret = Math.Max(ret, path[0, j]);
            for (int y = 1; y < lena; y++)
                for (int x = 1; x < lenb; x++)
                    ret = Math.Max(ret, path[y, x]);
            ret = Math.Max(ret, 0);
        }

        return ret;
    }
}
