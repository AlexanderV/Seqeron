namespace Seqeron.Genomics.Alignment;

// ksw_exts2_sse of minimap2 (ksw2_exts2_sse.c, SSE4.1 code path), lane by lane, on the same byte layout as the C `mem` (u, v, x, y, x2, donor, acceptor, s, sf = target, qr = reversed query).
public static partial class Ksw2Aligner
{
    /// <summary>
    /// <c>ksw_exts2_sse</c>: the spliced extension / global alignment of minimap2 (<c>-x splice</c>): affine gap cost <c>q + k e</c> for short gaps, a long gap (an intron) of cost <c>q2</c> plus the donor / acceptor score
    /// (<c>noncan</c> for a non-canonical site, the miniprot model with <see cref="Ksw2Option.SpliceComplex"/>, a bonus <paramref name="juncBonus"/> for a junction of the BED12 list <paramref name="junc"/>), no band; the
    /// introns are <c>N</c> operations of the CIGAR (those of at least <c>long_thres = (q2 - q) / e - 1 (+1)</c> bases).
    /// </summary>
    /// <param name="query">The query (symbols below <paramref name="m"/>; 0..3 = ACGT).</param>
    /// <param name="target">The target.</param>
    /// <param name="m">Alphabet size.</param>
    /// <param name="mat">The score matrix (<c>mat[0]</c> match, <c>mat[1]</c> mismatch, <c>mat[m*m-1]</c> wildcard).</param>
    /// <param name="q">Gap open penalty.</param>
    /// <param name="e">Gap extension penalty.</param>
    /// <param name="q2">Open penalty of an intron (must exceed <c>q + e</c>, otherwise the result is empty).</param>
    /// <param name="noncan">Penalty of a non-canonical splice site.</param>
    /// <param name="zdrop">Z-drop (negative: none).</param>
    /// <param name="juncBonus">Bonus of an annotated junction.</param>
    /// <param name="flag">The <see cref="Ksw2Option"/> flags (<see cref="Ksw2Option.SpliceForward"/> / <see cref="Ksw2Option.SpliceReverse"/>, <see cref="Ksw2Option.SpliceFlank"/>, <see cref="Ksw2Option.SpliceComplex"/>, <see cref="Ksw2Option.ReverseCigar"/> ...).</param>
    /// <param name="junc">Junction flags of the target positions (bit 0: forward donor, bit 1: forward acceptor, bit 2: reverse acceptor, bit 3: reverse donor), or null.</param>
    public static Ksw2Result ExtendS2(byte[] query, byte[] target, int m, sbyte[] mat, int q, int e, int q2, int noncan, int zdrop, int juncBonus, Ksw2Option flag = Ksw2Option.None, byte[]? junc = null)
    {
        CheckInputs(query, target, m, mat);
        CheckGapPenalties(q, e);
        int qlen = query.Length, tlen = target.Length;
        q = (sbyte)q;
        e = (sbyte)e;
        q2 = (sbyte)q2;
        noncan = (sbyte)noncan;
        juncBonus = (sbyte)juncBonus;
        var ez = new Ez();
        if (m <= 1 || qlen <= 0 || tlen <= 0 || q2 <= q + e)
            return ez.ToResult();
        if ((flag & Ksw2Option.SpliceForward) != 0 && (flag & Ksw2Option.SpliceReverse) != 0)
            throw new ArgumentException("The splice directions cannot both be set (assert in ksw_exts2_sse).", nameof(flag));
        int qe = q + e;
        bool withCigar = (flag & Ksw2Option.ScoreOnly) == 0;
        bool approxMax = (flag & Ksw2Option.ApproxMax) != 0;
        byte scMch = Sat(mat[0]);
        byte scMis = Sat(mat[1]);
        byte scN = mat[(m * m) - 1] == 0 ? Sat(-e) : Sat(mat[(m * m) - 1]);
        byte m1 = Sat(m - 1);
        byte qv = Sat(q), q2v = Sat(q2), qev = Sat(q + e);
        int tlen_ = (tlen + 15) / 16;
        int nCol = (((qlen < tlen ? qlen : tlen) + 15) / 16) + 1;
        int qlen_ = (qlen + 15) / 16;
        int maxSc = mat[0], minSc = mat[1];
        for (int t = 1; t < m * m; ++t)
        {
            maxSc = Math.Max(maxSc, mat[t]);
            minSc = Math.Min(minSc, mat[t]);
        }

        if (-minSc > 2 * (q + e))
            return ez.ToResult();
        int longThres = ((q2 - q) / e) - 1;
        if (q2 > q + e + (longThres * e))
            ++longThres;
        int longDiff = (longThres * e) - (q2 - q);

        var mem = new byte[((tlen_ * 9) + qlen_ + 1) * 16];
        int uO = 0, vO = tlen_ * 16, xO = 2 * tlen_ * 16, yO = 3 * tlen_ * 16, x2O = 4 * tlen_ * 16, donorO = 5 * tlen_ * 16, accO = 6 * tlen_ * 16, sO = 7 * tlen_ * 16, sfO = 8 * tlen_ * 16, qrO = sfO + (tlen_ * 16);
        Array.Fill(mem, Sat(-q - e), uO, tlen_ * 16 * 4);
        Array.Fill(mem, Sat(-q2), x2O, tlen_ * 16);
        int[]? H = null;
        if (!approxMax)
        {
            H = new int[tlen_ * 16];
            Array.Fill(H, NegInf);
        }

        byte[]? p = null;
        int[]? off = null, offEnd = null;
        if (withCigar)
        {
            RequireTraceback("ExtendS2", ((((long)(qlen + tlen - 1)) * nCol) + 1) * 16);
            p = new byte[((((long)(qlen + tlen - 1)) * nCol) + 1) * 16];
            off = new int[qlen + tlen - 1];
            offEnd = new int[qlen + tlen - 1];
        }

        for (int t = 0; t < qlen; ++t)
            mem[qrO + t] = query[qlen - 1 - t];
        Buffer.BlockCopy(target, 0, mem, sfO, tlen);

        // the donor and acceptor arrays (this assumes the 0/1/2/3 = A/C/G/T encoding)
        if ((flag & (Ksw2Option.SpliceForward | Ksw2Option.SpliceReverse)) != 0)
        {
            int[] sp0 = [8, 15, 21, 30];
            int[] sp = new int[4];
            if ((flag & Ksw2Option.SpliceComplex) != 0)
            {
                for (int t = 0; t < 4; ++t)
                    sp[t] = (int)((sp0[t] / 3.0) + .499);
            }
            else
            {
                sp[0] = (flag & Ksw2Option.SpliceFlank) != 0 ? noncan / 2 : 0;
                sp[1] = sp[2] = sp[3] = noncan;
            }

            Array.Fill(mem, Sat(-sp[3]), donorO, tlen_ * 16);
            Array.Fill(mem, Sat(-sp[3]), accO, tlen_ * 16);
            bool fwd = (flag & Ksw2Option.SpliceForward) != 0;
            if ((flag & Ksw2Option.ReverseCigar) == 0)
            {
                for (int t = 0; t < tlen - 4; ++t)
                {
                    int z = 3;
                    if (fwd)
                    {
                        if (target[t + 1] == 2 && target[t + 2] == 3)
                            z = target[t + 3] == 0 || target[t + 3] == 2 ? -1 : 0;
                        else if (target[t + 1] == 2 && target[t + 2] == 1)
                            z = 1;
                        else if (target[t + 1] == 0 && target[t + 2] == 3)
                            z = 2;
                    }
                    else
                    {
                        if (target[t + 1] == 1 && target[t + 2] == 3)
                            z = target[t + 3] == 0 || target[t + 3] == 2 ? -1 : 0;
                        else if (target[t + 1] == 2 && target[t + 2] == 3)
                            z = 2;
                    }

                    mem[donorO + t] = z < 0 ? (byte)0 : Sat(-sp[z]);
                }

                for (int t = 2; t < tlen; ++t)
                {
                    int z = 3;
                    if (fwd)
                    {
                        if (target[t - 1] == 0 && target[t] == 2)
                            z = target[t - 2] == 1 || target[t - 2] == 3 ? -1 : 0;
                        else if (target[t - 1] == 0 && target[t] == 1)
                            z = 2;
                    }
                    else
                    {
                        if (target[t - 1] == 0 && target[t] == 1)
                            z = target[t - 2] == 1 || target[t - 2] == 3 ? -1 : 0;
                        else if (target[t - 1] == 2 && target[t] == 1)
                            z = 1;
                        else if (target[t - 1] == 0 && target[t] == 3)
                            z = 2;
                    }

                    mem[accO + t] = z < 0 ? (byte)0 : Sat(-sp[z]);
                }
            }
            else
            {
                for (int t = 0; t < tlen - 4; ++t)
                {
                    int z = 3;
                    if (fwd)
                    {
                        if (target[t + 1] == 2 && target[t + 2] == 0)
                            z = target[t + 3] == 1 || target[t + 3] == 3 ? -1 : 0;
                        else if (target[t + 1] == 1 && target[t + 2] == 0)
                            z = 2;
                    }
                    else
                    {
                        if (target[t + 1] == 1 && target[t + 2] == 0)
                            z = target[t + 3] == 1 || target[t + 3] == 3 ? -1 : 0;
                        else if (target[t + 1] == 1 && target[t + 2] == 2)
                            z = 1;
                        else if (target[t + 1] == 3 && target[t + 2] == 0)
                            z = 2;
                    }

                    mem[donorO + t] = z < 0 ? (byte)0 : Sat(-sp[z]);
                }

                for (int t = 2; t < tlen; ++t)
                {
                    int z = 3;
                    if (fwd)
                    {
                        if (target[t - 1] == 3 && target[t] == 2)
                            z = target[t - 2] == 0 || target[t - 2] == 2 ? -1 : 0;
                        else if (target[t - 1] == 1 && target[t] == 2)
                            z = 1;
                        else if (target[t - 1] == 3 && target[t] == 0)
                            z = 2;
                    }
                    else
                    {
                        if (target[t - 1] == 3 && target[t] == 1)
                            z = target[t - 2] == 0 || target[t - 2] == 2 ? -1 : 0;
                        else if (target[t - 1] == 3 && target[t] == 2)
                            z = 2;
                    }

                    mem[accO + t] = z < 0 ? (byte)0 : Sat(-sp[z]);
                }
            }
        }

        if (junc is not null)
        {
            bool f1 = (flag & Ksw2Option.SpliceForward) != 0, f2 = (flag & Ksw2Option.SpliceReverse) != 0;
            if ((flag & Ksw2Option.ReverseCigar) == 0)
            {
                for (int t = 0; t < tlen - 1; ++t)
                    if ((f1 && (junc[t + 1] & 1) != 0) || (f2 && (junc[t + 1] & 8) != 0))
                        mem[donorO + t] = Sat(mem[donorO + t] + juncBonus);
                for (int t = 0; t < tlen; ++t)
                    if ((f1 && (junc[t] & 2) != 0) || (f2 && (junc[t] & 4) != 0))
                        mem[accO + t] = Sat(mem[accO + t] + juncBonus);
            }
            else
            {
                for (int t = 0; t < tlen - 1; ++t)
                    if ((f1 && (junc[t + 1] & 2) != 0) || (f2 && (junc[t + 1] & 4) != 0))
                        mem[donorO + t] = Sat(mem[donorO + t] + juncBonus);
                for (int t = 0; t < tlen; ++t)
                    if ((f1 && (junc[t] & 1) != 0) || (f2 && (junc[t] & 8) != 0))
                        mem[accO + t] = Sat(mem[accO + t] + juncBonus);
            }
        }

        int H0 = 0, lastH0t = 0, lastSt = -1, lastEn = -1;
        var xt1 = new byte[16];
        var x2t1 = new byte[16];
        var vt1 = new byte[16];
        for (int r = 0; r < qlen + tlen - 1; ++r)
        {
            AlignmentCancellation.ThrowIfRequested();
            int st = 0, en = tlen - 1;
            int qrr = qrO + (qlen - 1 - r);
            if (st < r - qlen + 1)
                st = r - qlen + 1;
            if (en > r)
                en = r;
            int st0 = st, en0 = en;
            st = st / 16 * 16;
            en = (((en + 16) / 16) * 16) - 1;
            int x1, x21, v1;
            int boundary = 0;
            if (r == 0)
                boundary = -q - e;
            else if (r < longThres)
                boundary = -e;
            else if (r == longThres)
                boundary = longDiff;
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
                    x21 = -q2;
                    v1 = -q - e;
                }
            }
            else
            {
                x1 = -q - e;
                x21 = -q2;
                v1 = boundary;
            }

            if (en >= r)
            {
                mem[yO + r] = Sat(-q - e);
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

            uint x1c = (byte)x1, x21c = (byte)x21, v1c = (byte)v1;
            int st_ = st / 16, en_ = en / 16;
            int prBase = (r * nCol) - st_;
            if (withCigar)
            {
                off![r] = st;
                offEnd![r] = en;
            }

            bool right = (flag & Ksw2Option.Right) != 0;
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
                    int a2a = S(Sat(a2 + S(mem[accO + bo + k])));
                    int d = 0;
                    if (!withCigar)
                    {
                        z = Math.Max(z, a);
                        z = Math.Max(z, b);
                        z = Math.Max(z, a2a);
                    }
                    else if (!right)
                    {
                        d = a > z ? 1 : 0;
                        z = Math.Max(z, a);
                        d = b > z ? 2 : d;
                        z = Math.Max(z, b);
                        d = a2a > z ? 3 : d;
                        z = Math.Max(z, a2a);
                    }
                    else
                    {
                        d = z > a ? 0 : 1;
                        z = Math.Max(z, a);
                        d = z > b ? d : 2;
                        z = Math.Max(z, b);
                        d = z > a2a ? d : 3;
                        z = Math.Max(z, a2a);
                    }

                    byte zb = Sat(z);
                    mem[uO + bo + k] = Sat(zb - vt1[k]);
                    mem[vO + bo + k] = Sat(zb - ut);
                    byte tmp = Sat(zb - qv);
                    a = S(Sat(a - tmp));
                    b = S(Sat(b - tmp));
                    a2 = S(Sat(a2 - Sat(zb - q2v)));
                    int donor = S(mem[donorO + bo + k]);
                    if (!withCigar)
                    {
                        mem[xO + bo + k] = Sat(Math.Max(a, 0) - qev);
                        mem[yO + bo + k] = Sat(Math.Max(b, 0) - qev);
                        mem[x2O + bo + k] = Sat(Math.Max(a2, donor) - q2v);
                        continue;
                    }

                    if (!right)
                    {
                        mem[xO + bo + k] = Sat((a > 0 ? a : 0) - qev);
                        if (a > 0)
                            d |= 0x08;
                        mem[yO + bo + k] = Sat((b > 0 ? b : 0) - qev);
                        if (b > 0)
                            d |= 0x10;
                        mem[x2O + bo + k] = Sat(Math.Max(a2, donor) - q2v);
                        if (a2 > donor)
                            d |= 0x20;
                    }
                    else
                    {
                        mem[xO + bo + k] = Sat((a >= 0 ? a : 0) - qev);
                        if (a >= 0)
                            d |= 0x08;
                        mem[yO + bo + k] = Sat((b >= 0 ? b : 0) - qev);
                        if (b >= 0)
                            d |= 0x10;
                        mem[x2O + bo + k] = Sat(Math.Max(donor, a2) - q2v);
                        if (donor <= a2)
                            d |= 0x20;
                    }

                    p![((long)(prBase + t) * 16) + k] = (byte)d;
                }
            }

            if (!approxMax)
            {
                int maxH, maxT;
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

                if (ApplyZdrop(ez, true, maxH, r, maxT, zdrop, 0))
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

                if ((flag & Ksw2Option.ApproxDrop) != 0 && ApplyZdrop(ez, true, H0, r, lastH0t, zdrop, 0))
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
                ez.Cigar = Backtrack(true, rev, p!, off!, offEnd, nCol * 16, tlen - 1, qlen - 1, longThres);
            else if (ez.MaxT >= 0 && ez.MaxQ >= 0)
                ez.Cigar = Backtrack(true, rev, p!, off!, offEnd, nCol * 16, ez.MaxT, ez.MaxQ, longThres);
        }

        return ez.ToResult();
    }
}
