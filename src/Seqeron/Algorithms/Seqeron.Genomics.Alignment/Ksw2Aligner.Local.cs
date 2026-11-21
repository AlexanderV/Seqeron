namespace Seqeron.Genomics.Alignment;

// ksw_ll_qinit / ksw_ll_i16 of minimap2 (ksw2_ll_sse.c): Farrar's striped Smith-Waterman with 16-bit lanes, lane by lane (the signed saturating adds, the unsigned saturating subtractions, the signed maxima, the lazy-F loop with its early
// exit and the choice of the end coordinates from the striped layout, including the zero-scored padding lanes past the end of the query), so that the score AND the end coordinates are those of the program.
public static partial class Ksw2Aligner
{
    /// <summary>The result of <see cref="LocalScore"/>: the best local alignment score and the 0-based query / target coordinates of its last base (-1 when the target is empty).</summary>
    public readonly record struct LocalScoreResult(int Score, int QueryEnd, int TargetEnd);

    private sealed class LocalQuery
    {
        internal int QLen, SLen, M;
        internal short[] Qp = [];
    }

    private static LocalQuery LocalQueryInit(int qlen, byte[] query, int m, sbyte[] mat)
    {
        const int p = 8;                       // size 2: 8 int16 lanes per __m128i
        int slen = (qlen + p - 1) / p;
        var q = new LocalQuery { QLen = qlen, SLen = slen, M = m, Qp = new short[slen * m * p] };
        int nlen = slen * p;
        int t = 0;
        for (int a = 0; a < m; ++a)
        {
            for (int i = 0; i < slen; ++i)
            {
                for (int k = i; k < nlen; k += slen)
                    q.Qp[t++] = (short)(k >= qlen ? 0 : mat[(a * m) + query[k]]);
            }
        }

        return q;
    }

    private static ushort AddsEpi16(ushort a, ushort b) => unchecked((ushort)Math.Clamp((short)a + (short)b, short.MinValue, short.MaxValue));

    private static ushort SubsEpu16(ushort a, ushort b) => a > b ? (ushort)(a - b) : (ushort)0;

    private static ushort MaxEpi16(ushort a, ushort b) => (short)a >= (short)b ? a : b;

    /// <summary>
    /// <c>ksw_ll_i16</c>: the local alignment score (affine gap cost <c>gapo + k * gape</c>) of the query against the target and the coordinates of its end, exactly as the 16-bit striped SSE2 function of ksw2 reports them
    /// (<c>ksw_ll_qinit</c> with size 2). The maximum is taken before the lazy-F correction as in the C code, the row of the best score is the LAST row that reaches it, and the query coordinate is the last cell of the
    /// striped array (segment-major) that holds the maximum.
    /// </summary>
    /// <param name="query">The query (symbols below <paramref name="m"/>).</param>
    /// <param name="target">The target.</param>
    /// <param name="m">Alphabet size.</param>
    /// <param name="mat">The <c>m * m</c> score matrix, indexed <c>mat[target * m + query]</c>.</param>
    /// <param name="gapOpen">Gap open penalty (a gap of length k costs <c>gapOpen + k * gapExtend</c>).</param>
    /// <param name="gapExtend">Gap extension penalty.</param>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveOptimization)]
    public static LocalScoreResult LocalScore(byte[] query, byte[] target, int m, sbyte[] mat, int gapOpen, int gapExtend)
    {
        CheckInputs(query, target, m, mat);
        var q = LocalQueryInit(query.Length, query, m, mat);
        int slen = q.SLen;
        int tlen = target.Length;
        var h0 = new ushort[slen * 8];
        var h1 = new ushort[slen * 8];
        var e = new ushort[slen * 8];
        var hmax = new ushort[slen * 8];
        ushort gapoe = unchecked((ushort)(gapOpen + gapExtend)), gape = unchecked((ushort)gapExtend);
        int gmax = 0, qe = -1, te = -1;
        var h = new ushort[8];
        var f = new ushort[8];
        var max = new ushort[8];
        for (int i = 0; i < tlen; ++i)
        {
            Array.Clear(f);
            Array.Clear(max);
            int sBase = target[i] * slen * 8;
            // h = _mm_slli_si128(H0[slen - 1], 2)
            for (int l = 1; l < 8; l++)
                h[l] = h0[((slen - 1) * 8) + l - 1];
            h[0] = 0;
            for (int j = 0; j < slen; ++j)
            {
                int vo = j * 8;
                for (int l = 0; l < 8; l++)
                {
                    ushort hv = AddsEpi16(h[l], unchecked((ushort)q.Qp[sBase + vo + l]));
                    ushort ee = e[vo + l];
                    hv = MaxEpi16(hv, ee);
                    hv = MaxEpi16(hv, f[l]);
                    max[l] = MaxEpi16(max[l], hv);
                    h1[vo + l] = hv;
                    hv = SubsEpu16(hv, gapoe);
                    ee = SubsEpu16(ee, gape);
                    ee = MaxEpi16(ee, hv);
                    e[vo + l] = ee;
                    f[l] = MaxEpi16(SubsEpu16(f[l], gape), hv);
                    h[l] = h0[vo + l];            // h = H0[j] for the next vector (read after the stores, as in C)
                }
            }

            bool done = false;
            for (int k = 0; k < 8 && !done; ++k)
            {
                for (int l = 7; l >= 1; l--)
                    f[l] = f[l - 1];
                f[0] = 0;
                for (int j = 0; j < slen; ++j)
                {
                    int vo = j * 8;
                    bool any = false;
                    for (int l = 0; l < 8; l++)
                    {
                        ushort hv = MaxEpi16(h1[vo + l], f[l]);
                        h1[vo + l] = hv;
                        hv = SubsEpu16(hv, gapoe);
                        f[l] = SubsEpu16(f[l], gape);
                        if ((short)f[l] > (short)hv)
                            any = true;
                    }

                    if (!any)
                    {
                        done = true;
                        break;
                    }
                }
            }

            // __max_8: horizontal signed maximum, then the lane 0 as an unsigned 16-bit value
            ushort mx = max[0];
            for (int l = 1; l < 8; l++)
                mx = MaxEpi16(mx, max[l]);
            int imax = mx;
            if (imax >= gmax)
            {
                gmax = imax;
                te = i;
                Array.Copy(h1, hmax, hmax.Length);
            }

            (h0, h1) = (h1, h0);
        }

        for (int i = 0; i < slen * 8; ++i)
        {
            if (hmax[i] == gmax)
                qe = (i / 8) + (i % 8 * slen);
        }

        return new LocalScoreResult(gmax, qe, te);
    }
}
