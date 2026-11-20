using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

// EMBOSS 6.6 needle (audit A2-28, B09 F122): nucleus/embaln.c embAlignPathCalcWithEndGapPenalties + embAlignWalkNWMatrixUsingCompass, in single precision as the program computes
// (float matrices, float gap penalties rounded by ajRoundFloat to 16 mantissa bits), with the same boundary initialisation, the same end-gap rules (-endweight / -endopen / -endextend),
// the same tie order in the traceback (LEFT / DOWN are followed only while the gap-extension identity holds within U_FEPS), so the rows are the rows of the program.
public static partial class SequenceAligner
{
    /// <summary>
    /// EMBOSS <c>needle</c>: global alignment (Needleman-Wunsch, Gotoh's three matrices) in the arithmetic and with the conventions of the program: a gap of length k costs
    /// <c>gapOpen + (k - 1) * gapExtend</c> (positive penalties, <c>-gapopen</c> / <c>-gapextend</c>, rounded as EMBOSS rounds them to 16 mantissa bits), the scores are single-precision floats, and
    /// the end gaps are free unless <paramref name="endWeight"/> (<c>-endweight</c>): then an end gap costs <paramref name="endOpen"/> + (k - 1) * <paramref name="endExtend"/> (<c>-endopen</c> / <c>-endextend</c>).
    /// The rows (ties included) and the score are those of the installed needle 6.6.0.
    /// </summary>
    /// <param name="sequence1">First sequence (<c>-asequence</c>).</param>
    /// <param name="sequence2">Second sequence (<c>-bsequence</c>).</param>
    /// <param name="matrix">Residue-pair scores (<c>-datafile</c>); default EDNAFULL for nucleotide sequences, BLOSUM62 (EBLOSUM62) otherwise.</param>
    /// <param name="gapOpen">Cost of a gap of length one (<c>-gapopen</c>, default 10).</param>
    /// <param name="gapExtend">Cost of each further gap position (<c>-gapextend</c>, default 0.5).</param>
    /// <param name="endWeight"><c>-endweight</c>: apply the end-gap penalties (default false: end gaps are free).</param>
    /// <param name="endOpen"><c>-endopen</c>: cost of an end gap of length one (default 10; used only with <paramref name="endWeight"/>).</param>
    /// <param name="endExtend"><c>-endextend</c>: cost of each further end-gap position (default 0.5; used only with <paramref name="endWeight"/>).</param>
    /// <returns>The global alignment (full rows, coordinates 0 .. length); <see cref="AlignmentResult.ExactScore"/> is the single-precision score of the program, <see cref="AlignmentResult.Score"/> its rounding.</returns>
    public static AlignmentResult NeedleAlign(string sequence1, string sequence2, SubstitutionMatrix? matrix = null, double gapOpen = 10, double gapExtend = 0.5, bool endWeight = false, double endOpen = 10, double endExtend = 0.5)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        if (gapOpen < 1 || gapOpen > 100)
            throw new ArgumentOutOfRangeException(nameof(gapOpen), "-gapopen is between 1.0 and 100.0 in EMBOSS needle.");
        if (gapExtend < 0 || gapExtend > 10)
            throw new ArgumentOutOfRangeException(nameof(gapExtend), "-gapextend is between 0.0 and 10.0 in EMBOSS needle.");
        if (endOpen < 1 || endOpen > 100)
            throw new ArgumentOutOfRangeException(nameof(endOpen), "-endopen is between 1.0 and 100.0 in EMBOSS needle.");
        if (endExtend < 0 || endExtend > 10)
            throw new ArgumentOutOfRangeException(nameof(endExtend), "-endextend is between 0.0 and 10.0 in EMBOSS needle.");
        if (sequence1.Length == 0 || sequence2.Length == 0)
            throw new ArgumentException("EMBOSS needle cannot align an empty sequence.", sequence1.Length == 0 ? nameof(sequence1) : nameof(sequence2));
        matrix ??= IsNucleotideText(sequence1.ToUpperInvariant()) && IsNucleotideText(sequence2.ToUpperInvariant()) ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62;
        var (score, r1, r2) = NeedleCore(sequence1, sequence2, matrix, RoundFloat((float)gapOpen, 8), RoundFloat((float)gapExtend, 8), (float)endOpen, (float)endExtend, endWeight);
        return new AlignmentResult(r1, r2, (int)Math.Round(score, MidpointRounding.AwayFromZero), AlignmentType.Global, 0, 0, sequence1.Length, sequence2.Length) { ExactScore = score };
    }

    /// <summary>ajRoundFloat (ajmath.c): keeps FLT_MANT_DIG - nbits (at least 8) mantissa bits, a remainder of exactly one half rounds down.</summary>
    internal static float RoundFloat(float a, uint nbits)
    {
        int bitsUsed = 24 - (int)nbits;
        if (bitsUsed < 8)
            bitsUsed = 8;
        double v = a;
        if (v == 0 || double.IsNaN(v) || double.IsInfinity(v))
            return a;
        int e = Math.ILogB(v);
        int i = e + 1;                         // frexp: v = x * 2^i with 0.5 <= |x| < 1
        double x = Math.ScaleB(v, -i);
        double y = Math.ScaleB(x, bitsUsed);
        double w = Math.Truncate(y);
        double z = y - w;
        if (z > 0.5)
            w += 1.0;
        if (z < -0.5)
            w -= 1.0;
        double b = Math.ScaleB(w, -bitsUsed);
        return (float)Math.ScaleB(b, i);
    }

    private const float NeedleFloatEpsilon = 1.192e-6F;     // U_FEPS

    private static bool FloatEq(float a, float b) => (b - NeedleFloatEpsilon) < a && a < (b + NeedleFloatEpsilon);

    // the traceback is a line-by-line port of embAlignPathCalcWithEndGapPenalties: its decision chain repeats the same move in several branches, which is kept so that the tie order stays visibly the program's
#pragma warning disable S1871
    private static (float Score, string Row1, string Row2) NeedleCore(string s1, string s2, SubstitutionMatrix matrix, float gapopen, float gapextend, float endgapopen, float endgapextend, bool endweight)
    {
        int lena = s1.Length;
        int lenb = s2.Length;
        if ((long)lena * lenb > 400_000_000)
            throw new ArgumentException("Sequences too big for needle (use StretcherAlign).");
        var table = new EmbossScoreTable(matrix);
        int[] ca = table.Codes(s1);
        int[] cb = table.Codes(s2);
        float Sub(int x, int y) => table.Score(x, y);

        long n = (long)lena * lenb;
        var m = new float[n];
        var ix = new float[n];
        var iy = new float[n];
        var compass = new int[n];

        if (!endweight)
        {
            endgapopen = 0;
            endgapextend = 0;
        }

        float match = Sub(ca[0], cb[0]);
        ix[0] = -endgapopen - gapopen;
        iy[0] = -endgapopen - gapopen;
        m[0] = match;
        long cursor = 0;
        long cursorp;
        float testog, testeg;
        for (int ypos = 1; ypos < lena; ++ypos)
        {
            match = Sub(ca[ypos], cb[0]);
            cursor = (long)ypos * lenb;
            cursorp = (long)(ypos - 1) * lenb;
            testog = m[cursorp] - gapopen;
            testeg = iy[cursorp] - gapextend;
            iy[cursor] = testog >= testeg ? testog : testeg;
            m[cursor] = match - (endgapopen + ((ypos - 1) * endgapextend));
            ix[cursor] = -endgapopen - (ypos * endgapextend) - gapopen;
        }

        ix[cursor] -= endgapopen;
        ix[cursor] += gapopen;
        cursor = 0;
        for (int xpos = 1; xpos < lenb; ++xpos)
        {
            match = Sub(ca[0], cb[xpos]);
            cursor = xpos;
            cursorp = xpos - 1;
            testog = m[cursorp] - gapopen;
            testeg = ix[cursorp] - gapextend;
            ix[cursor] = testog >= testeg ? testog : testeg;
            m[cursor] = match - (endgapopen + ((xpos - 1) * endgapextend));
            iy[cursor] = -endgapopen - (xpos * endgapextend) - gapopen;
        }

        iy[cursor] -= endgapopen;
        iy[cursor] += gapopen;

        int xp = 1;
        while (xp != lenb)
        {
            int yp = 1;
            int bconv = cb[xp];
            cursorp = xp - 1;
            cursor = xp++;
            while (yp < lena)
            {
                match = Sub(ca[yp++], bconv);
                cursor += lenb;
                float mp = m[cursorp];
                float ixp = ix[cursorp];
                float iyp = iy[cursorp];
                if (mp > ixp && mp > iyp)
                    m[cursor] = mp + match;
                else if (ixp > iyp)
                    m[cursor] = ixp + match;
                else
                    m[cursor] = iyp + match;

                if (xp == lenb)
                {
                    testog = m[++cursorp] - endgapopen;
                    testeg = iy[cursorp] - endgapextend;
                }
                else
                {
                    testog = m[++cursorp];
                    if (testog < ix[cursorp])
                        testog = ix[cursorp];
                    testog -= gapopen;
                    testeg = iy[cursorp] - gapextend;
                }

                iy[cursor] = testog > testeg ? testog : testeg;
                cursorp += lenb;
                if (yp == lena)
                {
                    testog = m[--cursorp] - endgapopen;
                    testeg = ix[cursorp] - endgapextend;
                }
                else
                {
                    testog = m[--cursorp];
                    if (testog < iy[cursorp])
                        testog = iy[cursorp];
                    testog -= gapopen;
                    testeg = ix[cursorp] - gapextend;
                }

                ix[cursor] = testog > testeg ? testog : testeg;
            }
        }

        // embAlignGetScoreNWMatrix with endweight = true (the program forces it): the optimum is in the last cell
        cursor = n - 1;
        float score;
        if (m[cursor] > ix[cursor] && m[cursor] > iy[cursor])
            score = m[cursor];
        else if (ix[cursor] > iy[cursor])
            score = ix[cursor];
        else
            score = iy[cursor];

        // traceback
        int xpos2 = lenb - 1;
        int ypos2 = lena - 1;
        int cp = 0;
        while (xpos2 >= 0 && ypos2 >= 0)
        {
            cursor = ((long)ypos2 * lenb) + xpos2;
            float mp = m[cursor];
            float extendHere = (ypos2 == 0 || ypos2 == lena - 1) ? endgapextend : gapextend;
            float extendHereB = (xpos2 == 0 || xpos2 == lenb - 1) ? endgapextend : gapextend;
            if (cp == 1 && FloatEq(extendHere, ix[cursor] - ix[cursor + 1]))
            {
                compass[cursor] = 1;
                xpos2--;
            }
            else if (cp == 2 && FloatEq(extendHereB, iy[cursor] - iy[cursor + lenb]))
            {
                compass[cursor] = 2;
                ypos2--;
            }
            else if (mp >= ix[cursor] && mp >= iy[cursor])
            {
                if (cp == 1 && FloatEq(mp, ix[cursor]))
                {
                    compass[cursor] = 1;
                    xpos2--;
                }
                else if (cp == 2 && FloatEq(mp, iy[cursor]))
                {
                    compass[cursor] = 2;
                    ypos2--;
                }
                else
                {
                    compass[cursor] = 0;
                    ypos2--;
                    xpos2--;
                }
            }
            else if (ix[cursor] >= iy[cursor] && xpos2 > -1)
            {
                compass[cursor] = 1;
                xpos2--;
            }
            else if (ypos2 > -1)
            {
                compass[cursor] = 2;
                ypos2--;
            }
            else
            {
                throw new InvalidOperationException("something is seriously wrong in the traceback algorithm");
            }

            cp = compass[cursor];
        }

        // embAlignWalkNWMatrixUsingCompass (the walk starts at the last cell: no trailing overhang)
        var mrow = new StringBuilder(lena + lenb);
        var nrow = new StringBuilder(lena + lenb);
        int xw = lenb - 1;
        int yw = lena - 1;
        // the program re-reads the start from the score routine: start2 = lenb - 1, start1 = lena - 1
        while (xw >= 0 && yw >= 0)
        {
            long cur = ((long)yw * lenb) + xw;
            if (compass[cur] == 0)
            {
                mrow.Append(s1[yw--]);
                nrow.Append(s2[xw--]);
            }
            else if (compass[cur] == 1)
            {
                mrow.Append('-');
                nrow.Append(s2[xw--]);
            }
            else
            {
                mrow.Append(s1[yw--]);
                nrow.Append('-');
            }
        }

        for (; xw >= 0; xw--)
        {
            nrow.Append(s2[xw]);
            mrow.Append('-');
        }

        for (; yw >= 0; yw--)
        {
            mrow.Append(s1[yw]);
            nrow.Append('-');
        }

        return (score, Reverse(mrow), Reverse(nrow));
    }
#pragma warning restore S1871

    private static string Reverse(StringBuilder sb)
    {
        var a = new char[sb.Length];
        for (int i = 0; i < a.Length; i++)
            a[i] = sb[a.Length - 1 - i];
        return new string(a);
    }
}
