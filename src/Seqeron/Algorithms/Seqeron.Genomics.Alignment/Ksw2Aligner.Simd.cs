namespace Seqeron.Genomics.Alignment;

// ksw_extz2_sse and ksw_extd2_sse of minimap2 (ksw2_extz2_sse.c, ksw2_extd2_sse.c, SSE4.1 code path), lane by lane. The buffers are one byte array laid out exactly like the C `mem` (u, v, x, y, [x2, y2,] s, sf = target, qr = reversed
// query) so that the 16-lane blocks that run past the band read and write the same neighbouring bytes as the C code does.
public static partial class Ksw2Aligner
{
    // 16-lane helpers over a byte buffer: the lane arrays are plain byte[16]
    private static byte Sat(int v) => unchecked((byte)v);

    private static sbyte S(byte b) => unchecked((sbyte)b);

    // slli_si128(x, 1) | carry: lane 0 = carry byte, lane k = x[k - 1]
    private static void ShiftIn(byte[] mem, int off, byte[] dst, uint carry)
    {
        dst[0] = (byte)(carry & 0xff);
        dst[1] = (byte)((carry >> 8) & 0xff);
        dst[2] = (byte)((carry >> 16) & 0xff);
        dst[3] = (byte)((carry >> 24) & 0xff);
        // lanes 1..3 of the carry are 0 for a byte value, 0xff for the sign extension of a negative int8 (_mm_cvtsi32_si128(int8)); they are OR-ed over the shifted lanes 1..3
        for (int k = 1; k < 16; k++)
            dst[k] = (byte)(mem[off + k - 1] | (k < 4 ? dst[k] : 0));
    }

    /// <summary>
    /// <c>ksw_extz2_sse</c>: extension / global alignment with affine gap cost <c>q + k e</c> by the anti-diagonal 8-bit difference recurrence, band <paramref name="w"/>, Z-drop, end bonus, ambiguous-base scoring
    /// (without <see cref="Ksw2Option.GenericScore"/>: <c>mat[0]</c> match, <c>mat[1]</c> mismatch, <c>mat[m*m-1]</c> wildcard for the last symbol <c>m - 1</c>, which is <c>-e</c> when 0) and CIGAR.
    /// </summary>
    public static Ksw2Result ExtendZ2(byte[] query, byte[] target, int m, sbyte[] mat, int q, int e, int w, int zdrop, int endBonus, Ksw2Option flag = Ksw2Option.None)
    {
        CheckInputs(query, target, m, mat);
        int qlen = query.Length;
        int tlen = target.Length;
        q = (sbyte)q;
        e = (sbyte)e;
        var ez = new Ez();
        if (m <= 0 || qlen <= 0 || tlen <= 0)
            return ez.ToResult();
        int qe = q + e;
        bool withCigar = (flag & Ksw2Option.ScoreOnly) == 0;
        bool approxMax = (flag & Ksw2Option.ApproxMax) != 0;
        byte scMch = Sat(mat[0]);
        byte scMis = Sat(mat[1]);
        byte scN = mat[(m * m) - 1] == 0 ? Sat(-e) : Sat(mat[(m * m) - 1]);
        byte m1 = Sat(m - 1);
        byte maxSc = Sat(mat[0] + ((q + e) * 2));
        byte qv = Sat(q);
        byte qe2v = Sat((q + e) * 2);
        if (w < 0)
            w = tlen > qlen ? tlen : qlen;
        int wl = w;
        int wr = w;
        int tlen_ = (tlen + 15) / 16;
        int nCol = qlen < tlen ? qlen : tlen;
        nCol = (((nCol < w + 1 ? nCol : w + 1) + 15) / 16) + 1;
        int qlen_ = (qlen + 15) / 16;
        int maxScore = mat[0];
        int minScore = mat[1];
        for (int t = 1; t < m * m; ++t)
        {
            maxScore = maxScore > mat[t] ? maxScore : mat[t];
            minScore = minScore < mat[t] ? minScore : mat[t];
        }

        if (-minScore > 2 * (q + e))
            return ez.ToResult();

        var mem = new byte[((tlen_ * 6) + qlen_ + 1) * 16];
        int uO = 0;
        int vO = tlen_ * 16;
        int xO = 2 * tlen_ * 16;
        int yO = 3 * tlen_ * 16;
        int sO = 4 * tlen_ * 16;
        int sfO = 5 * tlen_ * 16;
        int qrO = sfO + (tlen_ * 16);
        int[]? H = null;
        if (!approxMax)
        {
            H = new int[tlen_ * 16];
            Array.Fill(H, NegInf);
        }

        byte[]? p = null;
        int[]? off = null;
        int[]? offEnd = null;
        if (withCigar)
        {
            p = new byte[((((long)(qlen + tlen - 1)) * nCol) + 1) * 16];
            off = new int[qlen + tlen - 1];
            offEnd = new int[qlen + tlen - 1];
        }

        for (int t = 0; t < qlen; ++t)
            mem[qrO + t] = query[qlen - 1 - t];
        Buffer.BlockCopy(target, 0, mem, sfO, tlen);

        int H0 = 0;
        int lastH0t = 0;
        int lastSt = -1;
        int lastEn = -1;
        var xt1 = new byte[16];
        var vt1 = new byte[16];
        for (int r = 0; r < qlen + tlen - 1; ++r)
        {
            int st = 0;
            int en = tlen - 1;
            int qrr = qrO + (qlen - 1 - r);
            if (st < r - qlen + 1)
                st = r - qlen + 1;
            if (en > r)
                en = r;
            if (st < ((r - wr + 1) >> 1))
                st = (r - wr + 1) >> 1;
            if (en > ((r + wl) >> 1))
                en = (r + wl) >> 1;
            if (st > en)
            {
                ez.ZDropped = true;
                break;
            }

            int st0 = st;
            int en0 = en;
            st = st / 16 * 16;
            en = (((en + 16) / 16) * 16) - 1;
            int x1;
            int v1;
            if (st > 0)
            {
                if (st - 1 >= lastSt && st - 1 <= lastEn)
                {
                    x1 = S(mem[xO + st - 1]);
                    v1 = S(mem[vO + st - 1]);
                }
                else
                {
                    x1 = v1 = 0;
                }
            }
            else
            {
                x1 = 0;
                v1 = r != 0 ? S(qv) : 0;
            }

            if (en >= r)
            {
                mem[yO + r] = 0;
                mem[uO + r] = r != 0 ? qv : (byte)0;
            }

            // loop fission: the scores first
            if ((flag & Ksw2Option.GenericScore) == 0)
            {
                for (int t = st0; t <= en0; t += 16)
                {
                    for (int k = 0; k < 16; k++)
                    {
                        byte sq = mem[sfO + t + k];
                        byte stt = mem[qrr + t + k];
                        bool mask = sq == m1 || stt == m1;
                        byte tmp = sq == stt ? scMch : scMis;
                        mem[sO + t + k] = mask ? scN : tmp;
                    }
                }
            }
            else
            {
                for (int t = st0; t <= en0; ++t)
                    mem[sO + t] = unchecked((byte)mat[(mem[sfO + t] * m) + mem[qrr + t]]);
            }

            uint x1c = unchecked((uint)x1);      // _mm_cvtsi32_si128(int8): sign-extended into the low 32 bits
            uint v1c = unchecked((uint)v1);
            int st_ = st / 16;
            int en_ = en / 16;
            int prBase = (r * nCol) - st_;
            if (withCigar)
            {
                off![r] = st;
                offEnd![r] = en;
            }

            for (int t = st_; t <= en_; ++t)
            {
                int bo = t * 16;
                // block1
                ShiftIn(mem, xO + bo, xt1, x1c);
                x1c = mem[xO + bo + 15];
                ShiftIn(mem, vO + bo, vt1, v1c);
                v1c = mem[vO + bo + 15];
                for (int k = 0; k < 16; k++)
                {
                    byte z = Sat(mem[sO + bo + k] + qe2v);
                    byte a = Sat(xt1[k] + vt1[k]);
                    byte ut = mem[uO + bo + k];
                    byte b = Sat(mem[yO + bo + k] + ut);
                    int d = 0;
                    if (!withCigar)
                    {
                        z = (byte)Math.Max(S(z), S(a));
                    }
                    else if ((flag & Ksw2Option.Right) == 0)
                    {
                        d = S(a) > S(z) ? 1 : 0;
                        z = (byte)Math.Max(S(z), S(a));
                        d = S(b) > S(z) ? 2 : d;
                    }
                    else
                    {
                        d = S(z) > S(a) ? 0 : 1;
                        z = (byte)Math.Max(S(z), S(a));
                        d = S(z) > S(b) ? d : 2;
                    }

                    // block2
                    z = Math.Max(z, b);                       // _mm_max_epu8
                    z = Math.Min(z, maxSc);                   // _mm_min_epu8
                    mem[uO + bo + k] = Sat(z - vt1[k]);
                    mem[vO + bo + k] = Sat(z - ut);
                    z = Sat(z - qv);
                    a = Sat(a - z);
                    b = Sat(b - z);
                    if (!withCigar)
                    {
                        mem[xO + bo + k] = (byte)Math.Max((int)S(a), 0);
                        mem[yO + bo + k] = (byte)Math.Max((int)S(b), 0);
                        continue;
                    }

                    if ((flag & Ksw2Option.Right) == 0)
                    {
                        bool pa = S(a) > 0;
                        mem[xO + bo + k] = pa ? a : (byte)0;
                        if (pa)
                            d |= 0x08;
                        bool pb = S(b) > 0;
                        mem[yO + bo + k] = pb ? b : (byte)0;
                        if (pb)
                            d |= 0x10;
                    }
                    else
                    {
                        bool na = 0 > S(a);
                        mem[xO + bo + k] = na ? (byte)0 : a;
                        if (!na)
                            d |= 0x08;
                        bool nb = 0 > S(b);
                        mem[yO + bo + k] = nb ? (byte)0 : b;
                        if (!nb)
                            d |= 0x10;
                    }

                    p![((long)(prBase + t) * 16) + k] = (byte)d;
                }
            }

            if (!approxMax)
            {
                int maxH;
                int maxT;
                if (r > 0)
                {
                    maxH = H![en0] = en0 > 0 ? H[en0 - 1] + mem[uO + en0] - qe : H[en0] + mem[vO + en0] - qe;
                    maxT = en0;
                    int en1 = st0 + ((en0 - st0) / 4 * 4);
                    var maxHv = new int[4] { maxH, maxH, maxH, maxH };
                    var maxTv = new int[4] { maxT, maxT, maxT, maxT };
                    int t;
                    for (t = st0; t < en1; t += 4)
                    {
                        for (int k = 0; k < 4; k++)
                        {
                            int h1 = H[t + k] + mem[vO + t + k] - qe;
                            H[t + k] = h1;
                            if (h1 > maxHv[k])
                            {
                                maxHv[k] = h1;
                                maxTv[k] = t;
                            }
                        }
                    }

                    for (int i = 0; i < 4; ++i)
                    {
                        if (maxH < maxHv[i])
                        {
                            maxH = maxHv[i];
                            maxT = maxTv[i] + i;
                        }
                    }

                    for (; t < en0; ++t)
                    {
                        H[t] += mem[vO + t] - qe;
                        if (H[t] > maxH)
                        {
                            maxH = H[t];
                            maxT = t;
                        }
                    }
                }
                else
                {
                    H![0] = mem[vO] - qe - qe;
                    maxH = H[0];
                    maxT = 0;
                }

                if (en0 == tlen - 1 && H[en0] > ez.Mte)
                {
                    ez.Mte = H[en0];
                    ez.MteQ = r - en0;
                }

                if (r - st0 == qlen - 1 && H[st0] > ez.Mqe)
                {
                    ez.Mqe = H[st0];
                    ez.MqeT = st0;
                }

                if (ApplyZdrop(ez, true, maxH, r, maxT, zdrop, e))
                    break;
                if (r == qlen + tlen - 2 && en0 == tlen - 1)
                    ez.Score = H[tlen - 1];
            }
            else
            {
                if (r > 0)
                {
                    if (lastH0t >= st0 && lastH0t <= en0 && lastH0t + 1 >= st0 && lastH0t + 1 <= en0)
                    {
                        int d0 = mem[vO + lastH0t] - qe;
                        int d1 = mem[uO + lastH0t + 1] - qe;
                        if (d0 > d1)
                        {
                            H0 += d0;
                        }
                        else
                        {
                            H0 += d1;
                            ++lastH0t;
                        }
                    }
                    else if (lastH0t >= st0 && lastH0t <= en0)
                    {
                        H0 += mem[vO + lastH0t] - qe;
                    }
                    else
                    {
                        ++lastH0t;
                        H0 += mem[uO + lastH0t] - qe;
                    }

                    if ((flag & Ksw2Option.ApproxDrop) != 0 && ApplyZdrop(ez, true, H0, r, lastH0t, zdrop, e))
                        break;
                }
                else
                {
                    H0 = mem[vO] - qe - qe;
                    lastH0t = 0;
                }

                if (r == qlen + tlen - 2 && en0 == tlen - 1)
                    ez.Score = H0;
            }

            lastSt = st;
            lastEn = en;
        }

        if (withCigar)
        {
            bool rev = (flag & Ksw2Option.ReverseCigar) != 0;
            if (!ez.ZDropped && (flag & Ksw2Option.ExtensionOnly) == 0)
            {
                ez.Cigar = Backtrack(true, rev, p!, off!, offEnd, nCol * 16, tlen - 1, qlen - 1);
            }
            else if (!ez.ZDropped && (flag & Ksw2Option.ExtensionOnly) != 0 && ez.Mqe + endBonus > ez.Max)
            {
                ez.ReachEnd = true;
                ez.Cigar = Backtrack(true, rev, p!, off!, offEnd, nCol * 16, ez.MqeT, qlen - 1);
            }
            else if (ez.MaxT >= 0 && ez.MaxQ >= 0)
            {
                ez.Cigar = Backtrack(true, rev, p!, off!, offEnd, nCol * 16, ez.MaxT, ez.MaxQ);
            }
        }

        return ez.ToResult();
    }

    /// <summary>
    /// <c>ksw_extd2_sse</c>: the dual-affine gap cost <c>min(q + k e, q2 + k e2)</c> (the two pieces are swapped when <c>q2 + e2 &lt; q + e</c>), otherwise as <see cref="ExtendZ2"/>.
    /// </summary>
    public static Ksw2Result ExtendD2(byte[] query, byte[] target, int m, sbyte[] mat, int q, int e, int q2, int e2, int w, int zdrop, int endBonus, Ksw2Option flag = Ksw2Option.None)
    {
        CheckInputs(query, target, m, mat);
        int qlen = query.Length;
        int tlen = target.Length;
        q = (sbyte)q;
        e = (sbyte)e;
        q2 = (sbyte)q2;
        e2 = (sbyte)e2;
        var ez = new Ez();
        if (m <= 1 || qlen <= 0 || tlen <= 0)
            return ez.ToResult();
        if (q2 + e2 < q + e)
        {
            (q, q2) = (q2, q);
            (e, e2) = (e2, e);
        }

        int qe = q + e;
        bool withCigar = (flag & Ksw2Option.ScoreOnly) == 0;
        bool approxMax = (flag & Ksw2Option.ApproxMax) != 0;
        byte scMch = Sat(mat[0]);
        byte scMis = Sat(mat[1]);
        byte scN = mat[(m * m) - 1] == 0 ? Sat(-e2) : Sat(mat[(m * m) - 1]);
        byte m1 = Sat(m - 1);
        byte qv = Sat(q);
        byte q2v = Sat(q2);
        byte qev = Sat(q + e);
        byte qe2v = Sat(q2 + e2);
        if (w < 0)
            w = tlen > qlen ? tlen : qlen;
        int wl = w;
        int wr = w;
        int tlen_ = (tlen + 15) / 16;
        int nCol = qlen < tlen ? qlen : tlen;
        nCol = (((nCol < w + 1 ? nCol : w + 1) + 15) / 16) + 1;
        int qlen_ = (qlen + 15) / 16;
        int minScore = mat[1];
        for (int t = 1; t < m * m; ++t)
            minScore = minScore < mat[t] ? minScore : mat[t];
        if (-minScore > 2 * (q + e))
            return ez.ToResult();
        int longThres = e != e2 ? ((q2 - q) / (e - e2)) - 1 : 0;
        if (q2 + e2 + (longThres * e2) > q + e + (longThres * e))
            ++longThres;
        int longDiff = (longThres * (e - e2)) - (q2 - q) - e2;

        var mem = new byte[((tlen_ * 8) + qlen_ + 1) * 16];
        int uO = 0;
        int vO = tlen_ * 16;
        int xO = 2 * tlen_ * 16;
        int yO = 3 * tlen_ * 16;
        int x2O = 4 * tlen_ * 16;
        int y2O = 5 * tlen_ * 16;
        int sO = 6 * tlen_ * 16;
        int sfO = 7 * tlen_ * 16;
        int qrO = sfO + (tlen_ * 16);
        Array.Fill(mem, Sat(-q - e), uO, tlen_ * 16);
        Array.Fill(mem, Sat(-q - e), vO, tlen_ * 16);
        Array.Fill(mem, Sat(-q - e), xO, tlen_ * 16);
        Array.Fill(mem, Sat(-q - e), yO, tlen_ * 16);
        Array.Fill(mem, Sat(-q2 - e2), x2O, tlen_ * 16);
        Array.Fill(mem, Sat(-q2 - e2), y2O, tlen_ * 16);
        int[]? H = null;
        if (!approxMax)
        {
            H = new int[tlen_ * 16];
            Array.Fill(H, NegInf);
        }

        byte[]? p = null;
        int[]? off = null;
        int[]? offEnd = null;
        if (withCigar)
        {
            p = new byte[((((long)(qlen + tlen - 1)) * nCol) + 1) * 16];
            off = new int[qlen + tlen - 1];
            offEnd = new int[qlen + tlen - 1];
        }

        for (int t = 0; t < qlen; ++t)
            mem[qrO + t] = query[qlen - 1 - t];
        Buffer.BlockCopy(target, 0, mem, sfO, tlen);

        int H0 = 0;
        int lastH0t = 0;
        int lastSt = -1;
        int lastEn = -1;
        var xt1 = new byte[16];
        var x2t1 = new byte[16];
        var vt1 = new byte[16];
        for (int r = 0; r < qlen + tlen - 1; ++r)
        {
            int st = 0;
            int en = tlen - 1;
            int qrr = qrO + (qlen - 1 - r);
            if (st < r - qlen + 1)
                st = r - qlen + 1;
            if (en > r)
                en = r;
            if (st < ((r - wr + 1) >> 1))
                st = (r - wr + 1) >> 1;
            if (en > ((r + wl) >> 1))
                en = (r + wl) >> 1;
            if (st > en)
            {
                ez.ZDropped = true;
                break;
            }

            int st0 = st;
            int en0 = en;
            st = st / 16 * 16;
            en = (((en + 16) / 16) * 16) - 1;
            int x1;
            int x21;
            int v1;
            int boundary;
            if (r == 0)
                boundary = -q - e;
            else if (r < longThres)
                boundary = -e;
            else
                boundary = r == longThres ? longDiff : -e2;
            if (st > 0)
            {
                if (st - 1 >= lastSt && st - 1 <= lastEn)
                {
                    x1 = S(mem[xO + st - 1]);
                    x21 = S(mem[x2O + st - 1]);
                    v1 = S(mem[vO + st - 1]);
                }
                else
                {
                    x1 = -q - e;
                    x21 = -q2 - e2;
                    v1 = -q - e;
                }
            }
            else
            {
                x1 = -q - e;
                x21 = -q2 - e2;
                v1 = boundary;
            }

            if (en >= r)
            {
                mem[yO + r] = Sat(-q - e);
                mem[y2O + r] = Sat(-q2 - e2);
                mem[uO + r] = Sat(boundary);
            }

            if ((flag & Ksw2Option.GenericScore) == 0)
            {
                for (int t = st0; t <= en0; t += 16)
                {
                    for (int k = 0; k < 16; k++)
                    {
                        byte sq = mem[sfO + t + k];
                        byte stt = mem[qrr + t + k];
                        bool mask = sq == m1 || stt == m1;
                        byte tmp = sq == stt ? scMch : scMis;
                        mem[sO + t + k] = mask ? scN : tmp;
                    }
                }
            }
            else
            {
                for (int t = st0; t <= en0; ++t)
                    mem[sO + t] = unchecked((byte)mat[(mem[sfO + t] * m) + mem[qrr + t]]);
            }

            uint x1c = (byte)x1;                  // _mm_cvtsi32_si128((uint8_t) x1)
            uint x21c = (byte)x21;
            uint v1c = (byte)v1;
            int st_ = st / 16;
            int en_ = en / 16;
            int prBase = (r * nCol) - st_;
            if (withCigar)
            {
                off![r] = st;
                offEnd![r] = en;
            }

            for (int t = st_; t <= en_; ++t)
            {
                int bo = t * 16;
                ShiftIn(mem, xO + bo, xt1, x1c);
                x1c = mem[xO + bo + 15];
                ShiftIn(mem, vO + bo, vt1, v1c);
                v1c = mem[vO + bo + 15];
                ShiftIn(mem, x2O + bo, x2t1, x21c);
                x21c = mem[x2O + bo + 15];
                for (int k = 0; k < 16; k++)
                {
                    int z = S(mem[sO + bo + k]);
                    int a = S(Sat(xt1[k] + vt1[k]));
                    byte ut = mem[uO + bo + k];
                    int b = S(Sat(mem[yO + bo + k] + ut));
                    int a2 = S(Sat(x2t1[k] + vt1[k]));
                    int b2 = S(Sat(mem[y2O + bo + k] + ut));
                    int d = 0;
                    if (!withCigar)
                    {
                        z = Math.Max(z, a);
                        z = Math.Max(z, b);
                        z = Math.Max(z, a2);
                        z = Math.Max(z, b2);
                        z = Math.Min(z, S(scMch));
                    }
                    else if ((flag & Ksw2Option.Right) == 0)
                    {
                        d = a > z ? 1 : 0;
                        z = Math.Max(z, a);
                        d = b > z ? 2 : d;
                        z = Math.Max(z, b);
                        d = a2 > z ? 3 : d;
                        z = Math.Max(z, a2);
                        d = b2 > z ? 4 : d;
                        z = Math.Max(z, b2);
                        z = Math.Min(z, S(scMch));
                    }
                    else
                    {
                        d = z > a ? 0 : 1;
                        z = Math.Max(z, a);
                        d = z > b ? d : 2;
                        z = Math.Max(z, b);
                        d = z > a2 ? d : 3;
                        z = Math.Max(z, a2);
                        d = z > b2 ? d : 4;
                        z = Math.Max(z, b2);
                        z = Math.Min(z, S(scMch));
                    }

                    // block2
                    byte zb = Sat(z);
                    mem[uO + bo + k] = Sat(zb - vt1[k]);
                    mem[vO + bo + k] = Sat(zb - ut);
                    byte tmp = Sat(zb - qv);
                    a = S(Sat(a - tmp));
                    b = S(Sat(b - tmp));
                    tmp = Sat(zb - q2v);
                    a2 = S(Sat(a2 - tmp));
                    b2 = S(Sat(b2 - tmp));
                    if (!withCigar)
                    {
                        mem[xO + bo + k] = Sat(Math.Max(a, 0) - qev);
                        mem[yO + bo + k] = Sat(Math.Max(b, 0) - qev);
                        mem[x2O + bo + k] = Sat(Math.Max(a2, 0) - qe2v);
                        mem[y2O + bo + k] = Sat(Math.Max(b2, 0) - qe2v);
                        continue;
                    }

                    if ((flag & Ksw2Option.Right) == 0)
                    {
                        mem[xO + bo + k] = Sat((a > 0 ? a : 0) - qev);
                        if (a > 0)
                            d |= 0x08;
                        mem[yO + bo + k] = Sat((b > 0 ? b : 0) - qev);
                        if (b > 0)
                            d |= 0x10;
                        mem[x2O + bo + k] = Sat((a2 > 0 ? a2 : 0) - qe2v);
                        if (a2 > 0)
                            d |= 0x20;
                        mem[y2O + bo + k] = Sat((b2 > 0 ? b2 : 0) - qe2v);
                        if (b2 > 0)
                            d |= 0x40;
                    }
                    else
                    {
                        mem[xO + bo + k] = Sat(Math.Max(a, 0) - qev);
                        if (a >= 0)
                            d |= 0x08;
                        mem[yO + bo + k] = Sat(Math.Max(b, 0) - qev);
                        if (b >= 0)
                            d |= 0x10;
                        mem[x2O + bo + k] = Sat(Math.Max(a2, 0) - qe2v);
                        if (a2 >= 0)
                            d |= 0x20;
                        mem[y2O + bo + k] = Sat(Math.Max(b2, 0) - qe2v);
                        if (b2 >= 0)
                            d |= 0x40;
                    }

                    p![((long)(prBase + t) * 16) + k] = (byte)d;
                }
            }

            if (!approxMax)
            {
                int maxH;
                int maxT;
                if (r > 0)
                {
                    maxH = H![en0] = en0 > 0 ? H[en0 - 1] + S(mem[uO + en0]) : H[en0] + S(mem[vO + en0]);
                    maxT = en0;
                    int en1 = st0 + ((en0 - st0) / 4 * 4);
                    var maxHv = new int[4] { maxH, maxH, maxH, maxH };
                    var maxTv = new int[4] { maxT, maxT, maxT, maxT };
                    int t;
                    for (t = st0; t < en1; t += 4)
                    {
                        for (int k = 0; k < 4; k++)
                        {
                            int h1 = H[t + k] + S(mem[vO + t + k]);
                            H[t + k] = h1;
                            if (h1 > maxHv[k])
                            {
                                maxHv[k] = h1;
                                maxTv[k] = t;
                            }
                        }
                    }

                    for (int i = 0; i < 4; ++i)
                    {
                        if (maxH < maxHv[i])
                        {
                            maxH = maxHv[i];
                            maxT = maxTv[i] + i;
                        }
                    }

                    for (; t < en0; ++t)
                    {
                        H[t] += S(mem[vO + t]);
                        if (H[t] > maxH)
                        {
                            maxH = H[t];
                            maxT = t;
                        }
                    }
                }
                else
                {
                    H![0] = S(mem[vO]) - qe;
                    maxH = H[0];
                    maxT = 0;
                }

                if (en0 == tlen - 1 && H[en0] > ez.Mte)
                {
                    ez.Mte = H[en0];
                    ez.MteQ = r - en0;
                }

                if (r - st0 == qlen - 1 && H[st0] > ez.Mqe)
                {
                    ez.Mqe = H[st0];
                    ez.MqeT = st0;
                }

                if (ApplyZdrop(ez, true, maxH, r, maxT, zdrop, e2))
                    break;
                if (r == qlen + tlen - 2 && en0 == tlen - 1)
                    ez.Score = H[tlen - 1];
            }
            else
            {
                if (r > 0)
                {
                    if (lastH0t >= st0 && lastH0t <= en0 && lastH0t + 1 >= st0 && lastH0t + 1 <= en0)
                    {
                        int d0 = S(mem[vO + lastH0t]);
                        int d1 = S(mem[uO + lastH0t + 1]);
                        if (d0 > d1)
                        {
                            H0 += d0;
                        }
                        else
                        {
                            H0 += d1;
                            ++lastH0t;
                        }
                    }
                    else if (lastH0t >= st0 && lastH0t <= en0)
                    {
                        H0 += S(mem[vO + lastH0t]);
                    }
                    else
                    {
                        ++lastH0t;
                        H0 += S(mem[uO + lastH0t]);
                    }
                }
                else
                {
                    H0 = S(mem[vO]) - qe;
                    lastH0t = 0;
                }

                if ((flag & Ksw2Option.ApproxDrop) != 0 && ApplyZdrop(ez, true, H0, r, lastH0t, zdrop, e2))
                    break;
                if (r == qlen + tlen - 2 && en0 == tlen - 1)
                    ez.Score = H0;
            }

            lastSt = st;
            lastEn = en;
        }

        if (withCigar)
        {
            bool rev = (flag & Ksw2Option.ReverseCigar) != 0;
            if (!ez.ZDropped && (flag & Ksw2Option.ExtensionOnly) == 0)
            {
                ez.Cigar = Backtrack(true, rev, p!, off!, offEnd, nCol * 16, tlen - 1, qlen - 1);
            }
            else if (!ez.ZDropped && (flag & Ksw2Option.ExtensionOnly) != 0 && ez.Mqe + endBonus > ez.Max)
            {
                ez.ReachEnd = true;
                ez.Cigar = Backtrack(true, rev, p!, off!, offEnd, nCol * 16, ez.MqeT, qlen - 1);
            }
            else if (ez.MaxT >= 0 && ez.MaxQ >= 0)
            {
                ez.Cigar = Backtrack(true, rev, p!, off!, offEnd, nCol * 16, ez.MaxT, ez.MaxQ);
            }
        }

        return ez.ToResult();
    }
}
