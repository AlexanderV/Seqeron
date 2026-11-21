// The control flow of the ported C functions is kept literally (loop variables moved inside the body, empty search loops, the C compound conditions).
#pragma warning disable S127, S108, S1994, S3358, S1066, S1854

namespace Seqeron.Genomics.Alignment;

// align.c of minimap2 2.26 (mm_align_skeleton and everything it calls), ported function by function: the end-seed filters, the gap filling between the anchors with ksw2, the left and right extensions, the Z-drop and
// inversion tests, mm_fix_cigar, mm_update_extra, the =/X CIGAR, the inversion rescue and the re-estimation of the DP maxima. The kernels are Ksw2Aligner.ExtendZ2 / ExtendD2 / ExtendS2 / LocalScore.
public static partial class Minimap2Chaining
{
    private const ulong SeedLongJoin = 1UL << 40;
    private const ulong SeedIgnore = 1UL << 41;
    private const int CigarM = 0, CigarI = 1, CigarD = 2, CigarN = 3, CigarEq = 7, CigarX = 8;

    // ksw_extz_t
    private sealed class AlnEz
    {
        internal int Max, MaxQ = -1, MaxT = -1, Mqe = Ksw2Aligner.NegInf, MqeT = -1, Mte = Ksw2Aligner.NegInf, MteQ = -1, Score = Ksw2Aligner.NegInf;
        internal bool ZDropped, ReachEnd;
        internal List<uint> Cigar = [];

        // ksw_reset_extz
        internal void Reset()
        {
            MaxQ = MaxT = MqeT = MteQ = -1;
            Max = 0;
            Score = Mqe = Mte = Ksw2Aligner.NegInf;
            Cigar = [];
            ZDropped = false;
            ReachEnd = false;
        }

        internal void Load(Ksw2Result r)
        {
            Max = r.Max;
            ZDropped = r.ZDropped;
            MaxQ = r.MaxQ;
            MaxT = r.MaxT;
            Mqe = r.Mqe;
            MqeT = r.MqeT;
            Mte = r.Mte;
            MteQ = r.MteQ;
            Score = r.Score;
            ReachEnd = r.ReachEnd;
            Cigar = [.. r.Cigar.Select(c => ((uint)c.Length << 4) | (uint)"MIDN".IndexOf(c.Operation))];
        }
    }

    // ksw_gen_simple_mat(5, mat, a, b, sc_ambi) with its int8_t arguments
    private static sbyte[] GenSimpleMat(int a, int b, int scAmbi)
    {
        sbyte a8 = unchecked((sbyte)a), b8 = unchecked((sbyte)b), s8 = unchecked((sbyte)scAmbi);
        a8 = a8 < 0 ? unchecked((sbyte)-a8) : a8;
        b8 = b8 > 0 ? unchecked((sbyte)-b8) : b8;
        s8 = s8 > 0 ? unchecked((sbyte)-s8) : s8;
        const int m = 5;
        var mat = new sbyte[m * m];
        for (int i = 0; i < m - 1; ++i)
        {
            for (int j = 0; j < m - 1; ++j)
                mat[(i * m) + j] = i == j ? a8 : b8;
            mat[(i * m) + m - 1] = s8;
        }

        for (int j = 0; j < m; ++j)
            mat[((m - 1) * m) + j] = s8;
        return mat;
    }

    private static byte[] Rev(byte[] s, int off, int len)
    {
        var r = new byte[len];
        for (int i = 0; i < len; i++)
            r[i] = s[off + len - 1 - i];
        return r;
    }

    private static byte[] Slice(byte[] s, int off, int len)
    {
        var r = new byte[Math.Max(0, len)];
        Array.Copy(s, off, r, 0, r.Length);
        return r;
    }

    // a slice that reads outside the array as zeros (the C code reads the neighbouring bytes of the buffer)
    private static byte[] SliceZero(byte[] s, int off, int len)
    {
        var r = new byte[Math.Max(0, len)];
        for (int i = 0; i < r.Length; i++)
            if (off + i >= 0 && off + i < s.Length)
                r[i] = s[off + i];
        return r;
    }

    private static void UpdateMaxZdrop(int score, int i, int j, ref int max, ref int maxI, ref int maxJ, int e, ref int maxZdrop, int[,] pos)
    {
        if (score < max)
        {
            int li = i - maxI, lj = j - maxJ;
            int diff = li > lj ? li - lj : lj - li;
            int z = max - score - (diff * e);
            if (z > maxZdrop)
            {
                maxZdrop = z;
                pos[0, 0] = maxI;
                pos[0, 1] = i;
                pos[1, 0] = maxJ;
                pos[1, 1] = j;
            }
        }
        else
        {
            max = score;
            maxI = i;
            maxJ = j;
        }
    }

    // mm_test_zdrop: 0 no drop, 1 a Z-drop, 2 a Z-drop that might hide an inversion
    private static int TestZdrop(MinimapMapOptions o, int minDpMax, byte[] qseq, int qOff, byte[] tseq, IReadOnlyList<uint> cigar, sbyte[] mat)
    {
        var sc = o.Scoring;
        int score = 0, max = int.MinValue, maxI = -1, maxJ = -1, i = 0, j = 0, maxZdrop = 0;
        var pos = new int[2, 2] { { -1, -1 }, { -1, -1 } };
        foreach (uint c in cigar)
        {
            int op = (int)(c & 0xf), len = (int)(c >> 4);
            if (op == CigarM)
            {
                for (int l = 0; l < len; ++l)
                {
                    score += mat[(tseq[i + l] * 5) + qseq[qOff + j + l]];
                    UpdateMaxZdrop(score, i + l, j + l, ref max, ref maxI, ref maxJ, sc.GapExtend, ref maxZdrop, pos);
                }

                i += len;
                j += len;
            }
            else if (op is CigarI or CigarD or CigarN)
            {
                score -= sc.GapOpen + (sc.GapExtend * len);
                if (op == CigarI)
                    j += len;
                else
                    i += len;
                UpdateMaxZdrop(score, i, j, ref max, ref maxI, ref maxJ, sc.GapExtend, ref maxZdrop, pos);
            }
        }

        int qLen = pos[1, 1] - pos[1, 0], tLen = pos[0, 1] - pos[0, 0];
        if (!(o.Splice || o.ShortRead || o.Seed.ForwardOnly || o.Seed.ReverseOnly) && maxZdrop > sc.ZDropInverted && qLen < o.MaxGap && tLen < o.MaxGap)
        {
            var qseq2 = new byte[qLen];
            for (i = 0; i < qLen; ++i)
            {
                int c = qseq[qOff + pos[1, 1] - i - 1];
                qseq2[i] = c >= 4 ? (byte)4 : (byte)(3 - c);
            }

            var ll = Ksw2Aligner.LocalScore(qseq2, Slice(tseq, pos[0, 0], tLen), 5, mat, sc.GapOpen, sc.GapExtend);
            if (ll.Score >= o.MinChainScore * sc.Match && ll.Score >= minDpMax)
                return 2;
        }

        return maxZdrop > sc.ZDrop ? 1 : 0;
    }

    // mm_fix_cigar: left-align the indels, turn 5I6D7I into one I and one D, drop zero-length operations, remove a leading I / D
    private static void FixCigar(Reg r, byte[] qseq, int qBase, byte[] tseq, int tBase, out int qshift, out int tshift)
    {
        var cig = r.P!.Cigar;
        int toff = 0, qoff = 0;
        bool toShrink = false;
        qshift = tshift = 0;
        if (cig.Count <= 1)
            return;
        for (int k = 0; k < cig.Count; ++k)
        {
            uint op = cig[k] & 0xf, len = cig[k] >> 4;
            if (len == 0)
                toShrink = true;
            if (op == CigarM)
            {
                toff += (int)len;
                qoff += (int)len;
            }
            else if (op is CigarI or CigarD)
            {
                if (k > 0 && k < cig.Count - 1 && (cig[k - 1] & 0xf) == 0 && (cig[k + 1] & 0xf) == 0)
                {
                    int l, prevLen = (int)(cig[k - 1] >> 4);
                    if (op == CigarI)
                    {
                        for (l = 0; l < prevLen; ++l)
                            if (qseq[qBase + qoff - 1 - l] != qseq[qBase + qoff + (int)len - 1 - l])
                                break;
                    }
                    else
                    {
                        for (l = 0; l < prevLen; ++l)
                            if (tseq[tBase + toff - 1 - l] != tseq[tBase + toff + (int)len - 1 - l])
                                break;
                    }

                    if (l > 0)
                    {
                        cig[k - 1] -= (uint)l << 4;
                        cig[k + 1] += (uint)l << 4;
                        qoff -= l;
                        toff -= l;
                    }

                    if (l == prevLen)
                        toShrink = true;
                }

                if (op == CigarI)
                    qoff += (int)len;
                else
                    toff += (int)len;
            }
            else if (op == CigarN)
            {
                toff += (int)len;
            }
        }

        for (int k = 0; k < cig.Count - 2; ++k)
        {
            if ((cig[k] & 0xf) > 0 && (cig[k] & 0xf) + (cig[k + 1] & 0xf) == 3)
            {
                int l;
                var s = new uint[9];
                for (l = k; l < cig.Count; ++l)
                {
                    uint op = cig[l] & 0xf;
                    if (op is CigarI or CigarD || cig[l] >> 4 == 0)
                        s[op] += cig[l] >> 4;
                    else
                        break;
                }

                if (s[1] > 0 && s[2] > 0 && l - k > 2)
                {
                    cig[k] = (s[1] << 4) | CigarI;
                    cig[k + 1] = (s[2] << 4) | CigarD;
                    for (k += 2; k < l; ++k)
                        cig[k] &= 0xf;
                    toShrink = true;
                }

                k = l;
            }
        }

        if (toShrink)
        {
            var sq = new List<uint>(cig.Count);
            foreach (uint c in cig)
                if (c >> 4 != 0)
                    sq.Add(c);
            var merged = new List<uint>(sq.Count);
            for (int k = 0; k < sq.Count; ++k)
            {
                if (k == sq.Count - 1 || (sq[k] & 0xf) != (sq[k + 1] & 0xf))
                    merged.Add(sq[k]);
                else
                    sq[k + 1] += sq[k] >> 4 << 4;
            }

            cig = merged;
            r.P.Cigar = cig;
        }

        if ((cig[0] & 0xf) is CigarI or CigarD)
        {
            int l = (int)(cig[0] >> 4);
            if ((cig[0] & 0xf) == CigarI)
            {
                if (r.Rev)
                    r.Qe -= l;
                else
                    r.Qs += l;
                qshift = l;
            }
            else
            {
                r.Rs += l;
                tshift = l;
            }

            cig.RemoveAt(0);
        }
    }

    // mm_update_cigar_eqx: =/X instead of M
    private static void UpdateCigarEqx(Reg r, byte[] qseq, int qBase, byte[] tseq, int tBase)
    {
        if (r.P is null)
            return;
        int toff = 0, qoff = 0;
        var res = new List<uint>();
        foreach (uint c in r.P.Cigar)
        {
            int op = (int)(c & 0xf), len = (int)(c >> 4);
            if (op == CigarM)
            {
                while (len > 0)
                {
                    int l;
                    for (l = 0; l < len && qseq[qBase + qoff + l] == tseq[tBase + toff + l]; ++l)
                    {
                    }

                    if (l > 0)
                        res.Add(((uint)l << 4) | CigarEq);
                    len -= l;
                    toff += l;
                    qoff += l;
                    for (l = 0; l < len && qseq[qBase + qoff + l] != tseq[tBase + toff + l]; ++l)
                    {
                    }

                    if (l > 0)
                        res.Add(((uint)l << 4) | CigarX);
                    len -= l;
                    toff += l;
                    qoff += l;
                }

                continue;
            }

            if (op == CigarI)
                qoff += len;
            else if (op is CigarD or CigarN)
                toff += len;
            res.Add(c);
        }

        r.P.Cigar = res;
    }

    // mm_update_extra: left-align, count matches / block length / ambiguous bases and the maximal-scoring segment (dp_max)
    private static void UpdateExtra(Reg r, byte[] qseq, int qBase, byte[] tseq, int tBase, sbyte[] mat, int q, int e, bool isEqx, bool logGap)
    {
        var p = r.P;
        if (p is null)
            return;
        FixCigar(r, qseq, qBase, tseq, tBase, out int qshift, out int tshift);
        qBase += qshift;
        tBase += tshift;
        int toff = 0, qoff = 0;
        double s = 0.0, max = 0.0;
        r.Blen = r.Mlen = 0;
        foreach (uint c in p.Cigar)
        {
            int op = (int)(c & 0xf), len = (int)(c >> 4);
            if (op == CigarM)
            {
                int nAmbi = 0, nDiff = 0;
                for (int l = 0; l < len; ++l)
                {
                    int cq = qseq[qBase + qoff + l], ct = tseq[tBase + toff + l];
                    if (ct > 3 || cq > 3)
                        ++nAmbi;
                    else if (ct != cq)
                        ++nDiff;
                    s += mat[(ct * 5) + cq];
                    if (s < 0)
                        s = 0;
                    else
                        max = max > s ? max : s;
                }

                r.Blen += len - nAmbi;
                r.Mlen += len - (nAmbi + nDiff);
                p.NAmbi += nAmbi;
                toff += len;
                qoff += len;
            }
            else if (op == CigarI)
            {
                int nAmbi = 0;
                for (int l = 0; l < len; ++l)
                    if (qseq[qBase + qoff + l] > 3)
                        ++nAmbi;
                r.Blen += len - nAmbi;
                p.NAmbi += nAmbi;
                s -= logGap ? q + ((double)e * Log2((float)(1.0 + len))) : q + e;
                if (s < 0)
                    s = 0;
                qoff += len;
            }
            else if (op == CigarD)
            {
                int nAmbi = 0;
                for (int l = 0; l < len; ++l)
                    if (tseq[tBase + toff + l] > 3)
                        ++nAmbi;
                r.Blen += len - nAmbi;
                p.NAmbi += nAmbi;
                s -= logGap ? q + ((double)e * Log2((float)(1.0 + len))) : q + e;
                if (s < 0)
                    s = 0;
                toff += len;
            }
            else if (op == CigarN)
            {
                toff += len;
            }
        }

        p.DpMax = (int)(max + .499);
        if (isEqx)
            UpdateCigarEqx(r, qseq, qBase, tseq, tBase);
    }

    // mm_append_cigar
    private static void AppendCigar(Reg r, IReadOnlyList<uint> cigar)
    {
        if (cigar.Count == 0)
            return;
        r.P ??= new Extra();
        var c = r.P.Cigar;
        if (c.Count > 0 && (c[^1] & 0xf) == (cigar[0] & 0xf))
        {
            c[^1] += cigar[0] >> 4 << 4;
            for (int i = 1; i < cigar.Count; i++)
                c.Add(cigar[i]);
        }
        else
        {
            c.AddRange(cigar);
        }
    }

    // mm_align_pair
    private static void AlignPair(MinimapMapOptions o, byte[] qseq, byte[] tseq, byte[] junc, sbyte[] mat, int w, int endBonus, int zdrop, Ksw2Option flag, AlnEz ez)
    {
        var sc = o.Scoring;
        if (sc.MaxSmithWatermanMatrix > 0 && (long)tseq.Length * qseq.Length > sc.MaxSmithWatermanMatrix)
        {
            ez.Reset();
            ez.ZDropped = true;
        }
        else if (o.Splice)
        {
            var flagTmp = flag;
            if (!o.SpliceOld)
                flagTmp |= Ksw2Option.SpliceComplex;
            ez.Load(Ksw2Aligner.ExtendS2(qseq, tseq, 5, mat, sc.GapOpen, sc.GapExtend, sc.GapOpen2, sc.NonCanonical, zdrop, sc.JunctionBonus, flagTmp, junc));
        }
        else if (sc.GapOpen == sc.GapOpen2 && sc.GapExtend == sc.GapExtend2)
        {
            ez.Load(Ksw2Aligner.ExtendZ2(qseq, tseq, 5, mat, sc.GapOpen, sc.GapExtend, w, zdrop, endBonus, flag));
        }
        else
        {
            ez.Load(Ksw2Aligner.ExtendD2(qseq, tseq, 5, mat, sc.GapOpen, sc.GapExtend, sc.GapOpen2, sc.GapExtend2, w, zdrop, endBonus, flag));
        }
    }

    private static int GetHpLenBack(MmIndex mi, int rid, int x)
    {
        var seq = mi.Seq[rid];
        int off = x, i;
        int c = seq[off];
        for (i = off - 1; i >= 0; --i)
            if (seq[i] != c)
                break;
        return off - i;
    }

    // mm_adjust_minier
    private static void AdjustMinier(MmIndex mi, byte[][] qseq0, Mm128 a, out int r, out int q)
    {
        if (mi.Hpc)
        {
            var qseq = qseq0[a.X >> 63];
            q = (int)a.Y;
            int i, c;
            for (i = q - 1, c = qseq[q]; i > 0; --i)
                if (qseq[i] != c)
                    break;
            q = i + 1;
            c = GetHpLenBack(mi, (int)(a.X << 1 >> 33), (int)a.X);
            r = (int)a.X + 1 - c;
        }
        else
        {
            r = (int)a.X - (mi.K >> 1);
            q = (int)a.Y - (mi.K >> 1);
        }
    }

    private static int[]? CollectLongGaps(int as1, int cnt1, Mm128[] a, int minGap, out int nOut)
    {
        nOut = 0;
        int n = 0;
        for (int i = 1; i < cnt1; ++i)
        {
            int gap = ((int)a[as1 + i].Y - (int)a[as1 + i - 1].Y) - ((int)a[as1 + i].X - (int)a[as1 + i - 1].X);
            if (gap < -minGap || gap > minGap)
                ++n;
        }

        if (n <= 1)
            return null;
        var k = new int[n];
        n = 0;
        for (int i = 1; i < cnt1; ++i)
        {
            int gap = ((int)a[as1 + i].Y - (int)a[as1 + i - 1].Y) - ((int)a[as1 + i].X - (int)a[as1 + i - 1].X);
            if (gap < -minGap || gap > minGap)
                k[n++] = i;
        }

        nOut = n;
        return k;
    }

    // mm_filter_bad_seeds
    private static void FilterBadSeeds(int as1, int cnt1, Mm128[] a, int minGap, int diffThres, int maxExtLen, int maxExtCnt)
    {
        var kk = CollectLongGaps(as1, cnt1, a, minGap, out int n);
        if (kk is null)
            return;
        int maxSt = -1, maxEn = -1, max = 0;
        for (int k = 0; ; ++k)
        {
            int gap, l, nIns = 0, nDel = 0, qs, rs, maxDiff = 0, maxDiffL = -1;
            if (k == n || k >= maxEn)
            {
                if (maxEn > 0)
                    for (int i = kk[maxSt]; i < kk[maxEn]; ++i)
                        a[as1 + i].Y |= SeedIgnore;
                max = 0;
                maxSt = maxEn = -1;
                if (k == n)
                    break;
            }

            int ii = kk[k];
            gap = ((int)a[as1 + ii].Y - (int)a[as1 + ii - 1].Y) - (int)(a[as1 + ii].X - a[as1 + ii - 1].X);
            if (gap > 0)
                nIns += gap;
            else
                nDel += -gap;
            qs = (int)a[as1 + ii - 1].Y;
            rs = (int)a[as1 + ii - 1].X;
            for (l = k + 1; l < n && l <= k + maxExtCnt; ++l)
            {
                int j = kk[l], diff;
                if ((int)a[as1 + j].Y - qs > maxExtLen || (int)a[as1 + j].X - rs > maxExtLen)
                    break;
                gap = ((int)a[as1 + j].Y - (int)a[as1 + j - 1].Y) - (int)(a[as1 + j].X - a[as1 + j - 1].X);
                if (gap > 0)
                    nIns += gap;
                else
                    nDel += -gap;
                diff = nIns + nDel - Math.Abs(nIns - nDel);
                if (maxDiff < diff)
                {
                    maxDiff = diff;
                    maxDiffL = l;
                }
            }

            if (maxDiff > diffThres && maxDiff > max)
            {
                max = maxDiff;
                maxSt = k;
                maxEn = maxDiffL;
            }
        }
    }

    // mm_filter_bad_seeds_alt
    private static void FilterBadSeedsAlt(int as1, int cnt1, Mm128[] a, int minGap, int maxExt)
    {
        var kk = CollectLongGaps(as1, cnt1, a, minGap, out int n);
        if (kk is null)
            return;
        for (int k = 0; k < n;)
        {
            int i = kk[k], l;
            int gap1 = ((int)a[as1 + i].Y - (int)a[as1 + i - 1].Y) - ((int)a[as1 + i].X - (int)a[as1 + i - 1].X);
            int re1 = (int)a[as1 + i].X;
            int qe1 = (int)a[as1 + i].Y;
            gap1 = gap1 > 0 ? gap1 : -gap1;
            for (l = k + 1; l < n; ++l)
            {
                int j = kk[l], gap2, qSpanPre, rs2, qs2, m;
                if ((int)a[as1 + j].Y - qe1 > maxExt || (int)a[as1 + j].X - re1 > maxExt)
                    break;
                gap2 = ((int)a[as1 + j].Y - (int)a[as1 + j - 1].Y) - (int)(a[as1 + j].X - a[as1 + j - 1].X);
                qSpanPre = (int)(a[as1 + j - 1].Y >> 32 & 0xff);
                rs2 = (int)a[as1 + j - 1].X + qSpanPre;
                qs2 = (int)a[as1 + j - 1].Y + qSpanPre;
                m = rs2 - re1 < qs2 - qe1 ? rs2 - re1 : qs2 - qe1;
                gap2 = gap2 > 0 ? gap2 : -gap2;
                if (m > gap1 + gap2)
                    break;
                re1 = (int)a[as1 + j].X;
                qe1 = (int)a[as1 + j].Y;
                gap1 = gap2;
            }

            if (l > k + 1)
            {
                int end = kk[l - 1];
                for (int j = kk[k]; j < end; ++j)
                    a[as1 + j].Y |= SeedIgnore;
                a[as1 + end].Y |= SeedLongJoin;
            }

            k = l;
        }
    }

    // mm_fix_bad_ends
    private static void FixBadEnds(Reg r, Mm128[] a, int bw, int minMatch, out int asOut, out int cntOut)
    {
        asOut = r.As;
        cntOut = r.Cnt;
        if (r.Cnt < 3)
            return;
        int m, l;
        m = l = (int)(a[r.As].Y >> 32 & 0xff);
        for (int i = r.As + 1; i < r.As + r.Cnt - 1; ++i)
        {
            int qSpan = (int)(a[i].Y >> 32 & 0xff);
            if ((a[i].Y & SeedLongJoin) != 0)
                break;
            int lr = (int)a[i].X - (int)a[i - 1].X;
            int lq = (int)a[i].Y - (int)a[i - 1].Y;
            int min = lr < lq ? lr : lq;
            int max = lr > lq ? lr : lq;
            if (max - min > l >> 1)
                asOut = i;
            l += min;
            m += min < qSpan ? min : qSpan;
            if (l >= bw << 1 || (m >= minMatch && m >= bw) || m >= r.Mlen >> 1)
                break;
        }

        cntOut = r.As + r.Cnt - asOut;
        m = l = (int)(a[r.As + r.Cnt - 1].Y >> 32 & 0xff);
        for (int i = r.As + r.Cnt - 2; i > asOut; --i)
        {
            int qSpan = (int)(a[i + 1].Y >> 32 & 0xff);
            if ((a[i + 1].Y & SeedLongJoin) != 0)
                break;
            int lr = (int)a[i + 1].X - (int)a[i].X;
            int lq = (int)a[i + 1].Y - (int)a[i].Y;
            int min = lr < lq ? lr : lq;
            int max = lr > lq ? lr : lq;
            if (max - min > l >> 1)
                cntOut = i + 1 - asOut;
            l += min;
            m += min < qSpan ? min : qSpan;
            if (l >= bw << 1 || (m >= minMatch && m >= bw) || m >= r.Mlen >> 1)
                break;
        }
    }

    // mm_max_stretch
    private static void MaxStretch(Reg r, Mm128[] a, out int asOut, out int cntOut)
    {
        asOut = r.As;
        cntOut = r.Cnt;
        if (r.Cnt < 2)
            return;
        int maxScore = -1, maxI = -1, maxLen = 0;
        int score = (int)(a[r.As].Y >> 32 & 0xff), len = 1;
        int i;
        for (i = r.As + 1; i < r.As + r.Cnt; ++i)
        {
            int qSpan = (int)(a[i].Y >> 32 & 0xff);
            int lr = (int)a[i].X - (int)a[i - 1].X;
            int lq = (int)a[i].Y - (int)a[i - 1].Y;
            if (lq == lr)
            {
                score += lq < qSpan ? lq : qSpan;
                ++len;
            }
            else
            {
                if (score > maxScore)
                {
                    maxScore = score;
                    maxLen = len;
                    maxI = i - len;
                }

                score = qSpan;
                len = 1;
            }
        }

        if (score > maxScore)
        {
            maxScore = score;
            maxLen = len;
            maxI = i - len;
        }

        asOut = maxI;
        cntOut = maxLen;
    }

    // mm_seed_ext_score
    private static int SeedExtScore(MinimapMapOptions o, MmIndex mi, sbyte[] mat, int qlen, byte[][] qseq0, Mm128 a)
    {
        int qSpan = (int)(a.Y >> 32 & 0xff), extLen = o.AnchorExtensionLength;
        int rid = (int)(a.X << 1 >> 33);
        int re = (int)(uint)a.X + 1, rs = re - qSpan;
        int qe = (int)(uint)a.Y + 1, qs = qe - qSpan;
        rs = rs - extLen > 0 ? rs - extLen : 0;
        qs = qs - extLen > 0 ? qs - extLen : 0;
        re = re + extLen < mi.Length[rid] ? re + extLen : mi.Length[rid];
        qe = qe + extLen < qlen ? qe + extLen : qlen;
        byte[] qseq, tseq;
        if (o.Seed.QueryStrand)
        {
            qseq = Slice(qseq0[0], qs, qe - qs);
            tseq = mi.GetSeq2((a.X >> 63) != 0, rid, rs, re);
        }
        else
        {
            qseq = Slice(qseq0[a.X >> 63], qs, qe - qs);
            tseq = mi.GetSeq2(false, rid, rs, re);
        }

        return Ksw2Aligner.LocalScore(qseq, tseq, 5, mat, o.Scoring.GapOpen, o.Scoring.GapExtend).Score;
    }

    // mm_fix_bad_ends_splice
    private static void FixBadEndsSplice(MinimapMapOptions o, MmIndex mi, Reg r, sbyte[] mat, int qlen, byte[][] qseq0, Mm128[] a, out int as1, out int cnt1)
    {
        as1 = r.As;
        cnt1 = r.Cnt;
        if (r.Cnt < 3)
            return;
        double logGap = Math.Log((int)a[r.As + 1].X - (int)a[r.As].X);
        if ((a[r.As].Y >> 32 & 0xff) < logGap + o.AnchorExtensionShift)
        {
            int score = SeedExtScore(o, mi, mat, qlen, qseq0, a[r.As]);
            if ((double)score / mat[0] < logGap + o.AnchorExtensionShift)
            {
                ++as1;
                --cnt1;
            }
        }

        logGap = Math.Log((int)a[r.As + r.Cnt - 1].X - (int)a[r.As + r.Cnt - 2].X);
        if ((a[r.As + r.Cnt - 1].Y >> 32 & 0xff) < logGap + o.AnchorExtensionShift)
        {
            int score = SeedExtScore(o, mi, mat, qlen, qseq0, a[r.As + r.Cnt - 1]);
            if ((double)score / mat[0] < logGap + o.AnchorExtensionShift)
                --cnt1;
        }
    }

    // mm_split_reg
    private static Reg? SplitReg(Reg r, int n, int qlen, Mm128[] a, bool isQstrand)
    {
        if (n <= 0 || n >= r.Cnt)
            return null;
        var r2 = r.Clone();
        r2.Id = -1;
        r2.SamPri = false;
        r2.P = null;
        r2.SplitInv = false;
        r2.Cnt = r.Cnt - n;
        r2.Score = (int)(r.Score * ((float)r2.Cnt / r.Cnt) + .499);
        r2.As = r.As + n;
        if (r.Parent == r.Id)
            r2.Parent = ParentTmpPri;
        SetCoordinates(r2, qlen, a, isQstrand);
        r.Cnt -= r2.Cnt;
        r.Score -= r2.Score;
        SetCoordinates(r, qlen, a, isQstrand);
        r.Split |= 1;
        r2.Split |= 2;
        return r2;
    }

    // mm_align1: the base-level alignment of one chain; returns the second part when a Z-drop splits the chain
    private static Reg? Align1(MinimapMapOptions o, int minDpMax, MmIndex mi, int qlen, byte[][] qseq0, Reg r, int nA, Mm128[] a, AlnEz ez, bool spliceFor, bool spliceRev)
    {
        var sc = o.Scoring;
        bool isSr = o.ShortRead, isSplice = o.Splice;
        int rid = (int)(a[r.As].X << 1 >> 33);
        bool rev = (a[r.As].X >> 63) != 0;
        int as1, cnt1, rs0, re0, qs0, qe0, rs, re, qs, qe, rs1, qs1, re1, qe1;
        bool dropped = false;
        var extraFlag = Ksw2Option.None;
        if (r.Cnt == 0)
            return null;
        var mat = GenSimpleMat(sc.Match, sc.Mismatch, sc.Ambiguous);
        int bw = (int)((o.Bandwidth * 1.5) + 1.0);
        int bwLong = (int)((o.BandwidthLong * 1.5) + 1.0);
        if (bwLong < bw)
            bwLong = bw;

        if (isSr && !mi.Hpc)
        {
            MaxStretch(r, a, out as1, out cnt1);
            rs = (int)a[as1].X + 1 - (int)(a[as1].Y >> 32 & 0xff);
            qs = (int)a[as1].Y + 1 - (int)(a[as1].Y >> 32 & 0xff);
            re = (int)a[as1 + cnt1 - 1].X + 1;
            qe = (int)a[as1 + cnt1 - 1].Y + 1;
        }
        else
        {
            if (!o.NoEndFilter)
            {
                if (isSplice)
                    FixBadEndsSplice(o, mi, r, mat, qlen, qseq0, a, out as1, out cnt1);
                else
                    FixBadEnds(r, a, o.Bandwidth, o.MinChainScore * 2, out as1, out cnt1);
            }
            else
            {
                as1 = r.As;
                cnt1 = r.Cnt;
            }

            FilterBadSeeds(as1, cnt1, a, 10, 40, o.MaxGap >> 1, 10);
            FilterBadSeedsAlt(as1, cnt1, a, 30, o.MaxGap >> 1);
            AdjustMinier(mi, qseq0, a[as1], out rs, out qs);
            AdjustMinier(mi, qseq0, a[as1 + cnt1 - 1], out re, out qe);
        }

        if (isSplice)
        {
            if (spliceFor)
                extraFlag |= rev ? Ksw2Option.SpliceReverse : Ksw2Option.SpliceForward;
            if (spliceRev)
                extraFlag |= rev ? Ksw2Option.SpliceForward : Ksw2Option.SpliceReverse;
            if (o.SpliceFlank)
                extraFlag |= Ksw2Option.SpliceFlank;
        }

        int l;
        if (isSr)
        {
            qs0 = 0;
            qe0 = qlen;
            l = qs;
            l += l * sc.Match + sc.EndBonus > sc.GapOpen ? (l * sc.Match + sc.EndBonus - sc.GapOpen) / sc.GapExtend : 0;
            rs0 = rs - l > 0 ? rs - l : 0;
            l = qlen - qe;
            l += l * sc.Match + sc.EndBonus > sc.GapOpen ? (l * sc.Match + sc.EndBonus - sc.GapOpen) / sc.GapExtend : 0;
            re0 = re + l < mi.Length[rid] ? re + l : mi.Length[rid];
        }
        else
        {
            rs0 = (int)a[r.As].X + 1 - (int)(a[r.As].Y >> 32 & 0xff);
            qs0 = (int)a[r.As].Y + 1 - (int)(a[r.As].Y >> 32 & 0xff);
            if (rs0 < 0)
                rs0 = 0;
            rs1 = qs1 = 0;
            l = 0;
            for (int i = r.As - 1; i >= 0 && a[i].X >> 32 == a[r.As].X >> 32; --i)
            {
                int x = (int)a[i].X + 1 - (int)(a[i].Y >> 32 & 0xff);
                int y = (int)a[i].Y + 1 - (int)(a[i].Y >> 32 & 0xff);
                if (x < rs0 && y < qs0)
                {
                    if (++l > o.MinCount)
                    {
                        l = rs0 - x > qs0 - y ? rs0 - x : qs0 - y;
                        rs1 = rs0 - l;
                        qs1 = qs0 - l;
                        if (rs1 < 0)
                            rs1 = 0;
                        break;
                    }
                }
            }

            if (qs > 0 && rs > 0)
            {
                l = qs < o.MaxGap ? qs : o.MaxGap;
                qs1 = qs1 > qs - l ? qs1 : qs - l;
                qs0 = qs0 < qs1 ? qs0 : qs1;
                l += l * sc.Match > sc.GapOpen ? (l * sc.Match - sc.GapOpen) / sc.GapExtend : 0;
                l = l < o.MaxGap ? l : o.MaxGap;
                l = l < rs ? l : rs;
                rs1 = rs1 > rs - l ? rs1 : rs - l;
                rs0 = rs0 < rs1 ? rs0 : rs1;
                rs0 = rs0 < rs ? rs0 : rs;
            }
            else
            {
                rs0 = rs;
                qs0 = qs;
            }

            re0 = (int)a[r.As + r.Cnt - 1].X + 1;
            qe0 = (int)a[r.As + r.Cnt - 1].Y + 1;
            re1 = mi.Length[rid];
            qe1 = qlen;
            l = 0;
            for (int i = r.As + r.Cnt; i < nA && a[i].X >> 32 == a[r.As].X >> 32; ++i)
            {
                int x = (int)a[i].X + 1;
                int y = (int)a[i].Y + 1;
                if (x > re0 && y > qe0)
                {
                    if (++l > o.MinCount)
                    {
                        l = x - re0 > y - qe0 ? x - re0 : y - qe0;
                        re1 = re0 + l;
                        qe1 = qe0 + l;
                        break;
                    }
                }
            }

            if (qe < qlen && re < mi.Length[rid])
            {
                l = qlen - qe < o.MaxGap ? qlen - qe : o.MaxGap;
                qe1 = qe1 < qe + l ? qe1 : qe + l;
                qe0 = qe0 > qe1 ? qe0 : qe1;
                l += l * sc.Match > sc.GapOpen ? (l * sc.Match - sc.GapOpen) / sc.GapExtend : 0;
                l = l < o.MaxGap ? l : o.MaxGap;
                l = l < mi.Length[rid] - re ? l : mi.Length[rid] - re;
                re1 = re1 < re + l ? re1 : re + l;
                re0 = re0 > re1 ? re0 : re1;
            }
            else
            {
                re0 = re;
                qe0 = qe;
            }
        }

        if ((a[r.As].Y & SeedSelf) != 0)
        {
            int maxExt = r.Qs > r.Rs ? r.Qs - r.Rs : r.Rs - r.Qs;
            if (r.Rs - rs0 > maxExt)
                rs0 = r.Rs - maxExt;
            if (r.Qs - qs0 > maxExt)
                qs0 = r.Qs - maxExt;
            maxExt = r.Qe > r.Re ? r.Qe - r.Re : r.Re - r.Qe;
            if (re0 - r.Re > maxExt)
                re0 = r.Re + maxExt;
            if (qe0 - r.Qe > maxExt)
                qe0 = r.Qe + maxExt;
        }

        Reg? r2 = null;
        if (qs > 0 && rs > 0)
        {
            byte[] qseq, tseq, junc;
            if (o.Seed.QueryStrand)
            {
                qseq = Slice(qseq0[0], qs0, qs - qs0);
                tseq = mi.GetSeq2(rev, rid, rs0, rs);
            }
            else
            {
                qseq = Slice(qseq0[rev ? 1 : 0], qs0, qs - qs0);
                tseq = mi.GetSeq2(false, rid, rs0, rs);
            }

            junc = mi.BedJunction(rid, rs0, rs);
            qseq = Rev(qseq, 0, qs - qs0);
            tseq = Rev(tseq, 0, rs - rs0);
            junc = Rev(junc, 0, rs - rs0);
            AlignPair(o, qseq, tseq, junc, mat, bw, sc.EndBonus, r.SplitInv ? sc.ZDropInverted : sc.ZDrop, extraFlag | Ksw2Option.ExtensionOnly | Ksw2Option.Right | Ksw2Option.ReverseCigar, ez);
            if (ez.Cigar.Count > 0)
            {
                AppendCigar(r, ez.Cigar);
                r.P!.DpScore += ez.Max;
            }

            rs1 = rs - (ez.ReachEnd ? ez.MqeT + 1 : ez.MaxT + 1);
            qs1 = qs - (ez.ReachEnd ? qs - qs0 : ez.MaxQ + 1);
        }
        else
        {
            rs1 = rs;
            qs1 = qs;
        }

        re1 = rs;
        qe1 = qs;

        for (int i = isSr ? cnt1 - 1 : 1; i < cnt1; ++i)
        {
            if ((a[as1 + i].Y & (SeedIgnore | SeedTandem)) != 0 && i != cnt1 - 1)
                continue;
            if (isSr && !mi.Hpc)
            {
                re = (int)a[as1 + i].X + 1;
                qe = (int)a[as1 + i].Y + 1;
            }
            else
            {
                AdjustMinier(mi, qseq0, a[as1 + i], out re, out qe);
            }

            re1 = re;
            qe1 = qe;
            if (i == cnt1 - 1 || (a[as1 + i].Y & SeedLongJoin) != 0 || (qe - qs >= o.MinKswLength && re - rs >= o.MinKswLength))
            {
                int bw1 = bwLong, zdropCode;
                if ((a[as1 + i].Y & SeedLongJoin) != 0)
                    bw1 = qe - qs > re - rs ? qe - qs : re - rs;
                byte[] qseq, tseq;
                if (o.Seed.QueryStrand)
                {
                    qseq = Slice(qseq0[0], qs, qe - qs);
                    tseq = mi.GetSeq2(rev, rid, rs, re);
                }
                else
                {
                    qseq = Slice(qseq0[rev ? 1 : 0], qs, qe - qs);
                    tseq = mi.GetSeq2(false, rid, rs, re);
                }

                byte[] junc = mi.BedJunction(rid, rs, re);
                if (isSr)
                {
                    ez.Reset();
                    ez.Score = 0;
                    for (int j = 0; j < qe - qs; ++j)
                    {
                        if (qseq[j] >= 4 || tseq[j] >= 4)
                            ez.Score += sc.GapExtend2;
                        else
                            ez.Score += qseq[j] == tseq[j] ? sc.Match : -sc.Mismatch;
                    }

                    ez.Cigar = [((uint)(qe - qs) << 4) | CigarM];
                }
                else
                {
                    AlignPair(o, qseq, tseq, junc, mat, bw1, -1, sc.ZDrop, extraFlag | Ksw2Option.ApproxMax, ez);
                }

                if ((zdropCode = TestZdrop(o, minDpMax, qseq, 0, tseq, ez.Cigar, mat)) != 0)
                    AlignPair(o, qseq, tseq, junc, mat, bw1, -1, zdropCode == 2 ? sc.ZDropInverted : sc.ZDrop, extraFlag, ez);
                if (ez.Cigar.Count > 0)
                    AppendCigar(r, ez.Cigar);
                if (ez.ZDropped)
                {
                    r.P ??= new Extra();
                    int j;
                    for (j = i - 1; j >= 0; --j)
                        if ((int)a[as1 + j].X <= rs + ez.MaxT)
                            break;
                    dropped = true;
                    if (j < 0)
                        j = 0;
                    r.P.DpScore += ez.Max;
                    re1 = rs + (ez.MaxT + 1);
                    qe1 = qs + (ez.MaxQ + 1);
                    if (cnt1 - (j + 1) >= o.MinCount)
                    {
                        r2 = SplitReg(r, as1 + j + 1 - r.As, qlen, a, o.Seed.QueryStrand);
                        if (r2 is not null && zdropCode == 2)
                            r2.SplitInv = true;
                    }

                    break;
                }

                r.P!.DpScore += ez.Score;
                rs = re;
                qs = qe;
            }
        }

        if (!dropped && qe < qe0 && re < re0)
        {
            byte[] qseq, tseq;
            if (o.Seed.QueryStrand)
            {
                qseq = Slice(qseq0[0], qe, qe0 - qe);
                tseq = mi.GetSeq2(rev, rid, re, re0);
            }
            else
            {
                qseq = Slice(qseq0[rev ? 1 : 0], qe, qe0 - qe);
                tseq = mi.GetSeq2(false, rid, re, re0);
            }

            byte[] junc = mi.BedJunction(rid, re, re0);
            AlignPair(o, qseq, tseq, junc, mat, bw, sc.EndBonus, sc.ZDrop, extraFlag | Ksw2Option.ExtensionOnly, ez);
            if (ez.Cigar.Count > 0)
            {
                AppendCigar(r, ez.Cigar);
                r.P!.DpScore += ez.Max;
            }

            re1 = re + (ez.ReachEnd ? ez.MqeT + 1 : ez.MaxT + 1);
            qe1 = qe + (ez.ReachEnd ? qe0 - qe : ez.MaxQ + 1);
        }

        r.Rs = rs1;
        r.Re = re1;
        if (!rev || o.Seed.QueryStrand)
        {
            r.Qs = qs1;
            r.Qe = qe1;
        }
        else
        {
            r.Qs = qlen - qe1;
            r.Qe = qlen - qs1;
        }

        if (r.P is not null)
        {
            byte[] tseq = mi.GetSeq2(o.Seed.QueryStrand && r.Rev, rid, rs1, re1);
            var qarr = o.Seed.QueryStrand ? qseq0[0] : qseq0[r.Rev ? 1 : 0];
            UpdateExtra(r, qarr, qs1, tseq, 0, mat, sc.GapOpen, sc.GapExtend, o.Eqx, !isSr);
            if (rev && r.P.TransStrand != 0)
                r.P.TransStrand ^= 3;
        }

        return r2;
    }

    // mm_align1_inv: rescue an inversion between two hits of one split chain
    private static Reg? Align1Inv(MinimapMapOptions o, int minDpMax, MmIndex mi, int qlen, byte[][] qseq0, Reg r1, Reg r2, AlnEz ez)
    {
        var sc = o.Scoring;
        if ((r1.Split & 1) == 0 || (r2.Split & 2) == 0)
            return null;
        if (r1.Id != r1.Parent && r1.Parent != ParentTmpPri)
            return null;
        if (r2.Id != r2.Parent && r2.Parent != ParentTmpPri)
            return null;
        if (r1.Rid != r2.Rid || r1.Rev != r2.Rev)
            return null;
        int ql = r1.Rev ? r1.Qs - r2.Qe : r2.Qs - r1.Qe;
        int tl = r2.Rs - r1.Re;
        if (ql < o.MinChainScore || ql > o.MaxGap)
            return null;
        if (tl < o.MinChainScore || tl > o.MaxGap)
            return null;
        var mat = GenSimpleMat(sc.Match, sc.Mismatch, sc.Ambiguous);
        var tseq0 = mi.GetSeq(r1.Rid, r1.Re, r2.Rs) ?? [];
        // the program points into the whole query array: the end coordinate of the striped SW can lie in the zero-scored padding past the query (q_off < 0), then the alignment starts before the region
        var qarr = r1.Rev ? qseq0[0] : qseq0[1];
        int qbase = r1.Rev ? r2.Qe : qlen - r2.Qs;
        var qrev = Rev(qarr, qbase, ql);
        var trev = Rev(tseq0, 0, tl);
        var ll = Ksw2Aligner.LocalScore(qrev, trev, 5, mat, sc.GapOpen, sc.GapExtend);
        if (ll.Score < minDpMax)
            return null;
        int qOff = ql - (ll.QueryEnd + 1), tOff = tl - (ll.TargetEnd + 1);
        var qa = SliceZero(qarr, qbase + qOff, ql - qOff);
        var ta = Slice(tseq0, tOff, tl - tOff);
        AlignPair(o, qa, ta, [], mat, (int)(o.Bandwidth * 1.5), -1, sc.ZDrop, Ksw2Option.ExtensionOnly, ez);
        if (ez.Cigar.Count == 0)
            return null;
        var rInv = new Reg();
        AppendCigar(rInv, ez.Cigar);
        rInv.P!.DpScore = ez.Max;
        rInv.Id = -1;
        rInv.Parent = ParentUnset;
        rInv.Inv = true;
        rInv.Rev = !r1.Rev;
        rInv.Rid = r1.Rid;
        rInv.Div = -1.0f;
        if (!rInv.Rev)
        {
            rInv.Qs = r2.Qe + qOff;
            rInv.Qe = rInv.Qs + ez.MaxQ + 1;
        }
        else
        {
            rInv.Qe = r2.Qs - qOff;
            rInv.Qs = rInv.Qe - (ez.MaxQ + 1);
        }

        rInv.Rs = r1.Re + tOff;
        rInv.Re = rInv.Rs + ez.MaxT + 1;
        UpdateExtra(rInv, qarr, qbase + qOff, tseq0, tOff, mat, sc.GapOpen, sc.GapExtend, o.Eqx, !o.ShortRead);
        return rInv;
    }

    // mm_event_identity
    private static double EventIdentity(Reg r)
    {
        if (r.P is null)
            return -1.0;
        CountGaps(r, out int nGap, out int nGapo);
        return (double)r.Mlen / (r.Blen + r.P.NAmbi - nGap + nGapo);
    }

    private static void CountGaps(Reg r, out int nGap, out int nGapo)
    {
        nGap = nGapo = 0;
        if (r.P is null)
            return;
        foreach (uint c in r.P.Cigar)
        {
            int op = (int)(c & 0xf), len = (int)(c >> 4);
            if (op is CigarI or CigarD)
            {
                ++nGapo;
                nGap += len;
            }
        }
    }

    // mm_recal_max_dp
    private static int RecalMaxDp(Reg r, double b2, int matchSc)
    {
        if (r.P is null)
            return -1;
        int nGap = 0, nGapo = 0;
        double gapCost = 0.0;
        foreach (uint c in r.P.Cigar)
        {
            int op = (int)(c & 0xf), len = (int)(c >> 4);
            if (op is CigarI or CigarD)
            {
                gapCost += b2 + (double)Log2((float)(1.0 + len));
                ++nGapo;
                nGap += len;
            }
        }

        int nMis = r.Blen + r.P.NAmbi - r.Mlen - nGap;
        return (int)((matchSc * (r.Mlen - (b2 * nMis) - gapCost)) + .499);
    }

    // mm_update_dp_max
    private static void UpdateDpMax(int qlen, List<Reg> regs, float frac, int a, int b)
    {
        int max = -1, max2 = -1, maxI = -1;
        if (regs.Count < 2)
            return;
        for (int i = 0; i < regs.Count; ++i)
        {
            var r = regs[i];
            if (r.P is null)
                continue;
            if (r.P.DpMax > max)
            {
                max2 = max;
                max = r.P.DpMax;
                maxI = i;
            }
            else if (r.P.DpMax > max2)
            {
                max2 = r.P.DpMax;
            }
        }

        if (maxI < 0 || max < 0 || max2 < 0)
            return;
        if (regs[maxI].Qe - regs[maxI].Qs < (double)qlen * frac)
            return;
        if (max2 < (double)max * frac)
            return;
        double div = 1.0 - EventIdentity(regs[maxI]);
        if (div < 0.02)
            div = 0.02;
        double b2 = 0.5 / div;
        if (b2 * a < b)
            b2 = (double)a / b;
        foreach (var r in regs)
        {
            if (r.P is null)
                continue;
            r.P.DpMax = RecalMaxDp(r, b2, a);
            if (r.P.DpMax < 0)
                r.P.DpMax = 0;
        }
    }

    // mm_squeeze_a: move the anchors of the hits to the front of the array (in the order of their offsets) and return the new length
    private static int SqueezeA(List<Reg> regs, Mm128[] a)
    {
        int asPos = 0;
        foreach (var r in regs.OrderBy(x => x.As))
        {
            if (r.As != asPos)
            {
                Array.Copy(a, r.As, a, asPos, r.Cnt);
                r.As = asPos;
            }

            asPos += r.Cnt;
        }

        return asPos;
    }

    // mm_filter_regs
    private static void FilterRegs(MinimapMapOptions o, int minDpMax, int qlen, List<Reg> regs)
    {
        int k = 0;
        for (int i = 0; i < regs.Count; ++i)
        {
            var r = regs[i];
            bool flt = false;
            if (!r.Inv && !r.SegSplit && r.Cnt < o.MinCount)
                flt = true;
            if (r.P is not null)
            {
                if (r.Mlen < o.MinChainScore)
                    flt = true;
                else if (r.P.DpMax < minDpMax)
                    flt = true;
                else if (r.Qs > qlen * o.MaxClipRatio && qlen - r.Qe > qlen * o.MaxClipRatio)
                    flt = true;
            }

            if (!flt)
                regs[k++] = r;
        }

        regs.RemoveRange(k, regs.Count - k);
    }

    private static int MinDpMaxOf(MinimapMapOptions o) => o.Scoring.MinDpMax ?? 80;

    // mm_align_skeleton: align every hit; returns the new list of hits (the number can change: split chains, inversion hits, filtered hits)
    private static List<Reg> AlignSkeleton(MinimapMapOptions o, MmIndex mi, int qlen, string qstr, List<Reg> regs, Mm128[] a)
    {
        int minDpMax = MinDpMaxOf(o);
        var qseq0 = new byte[2][];
        qseq0[0] = new byte[qlen];
        qseq0[1] = new byte[qlen];
        for (int i = 0; i < qlen; ++i)
        {
            qseq0[0][i] = (byte)Nt4(qstr[i]);
            qseq0[1][qlen - 1 - i] = qseq0[0][i] < 4 ? (byte)(3 - qseq0[0][i]) : (byte)4;
        }

        int nA = SqueezeA(regs, a);
        var ez = new AlnEz();
        for (int i = 0; i < regs.Count; ++i)
        {
            Reg? r2;
            if (o.Splice && o.SpliceForward && o.SpliceReverse)
            {
                var s0 = regs[i].Clone();
                var s1 = regs[i].Clone();
                var s20 = Align1(o, minDpMax, mi, qlen, qseq0, s0, nA, a, ez, true, false);
                var s21 = Align1(o, minDpMax, mi, qlen, qseq0, s1, nA, a, ez, false, true);
                int which, transStrand;
                int d0 = s0.P?.DpScore ?? 0, d1 = s1.P?.DpScore ?? 0;
                if (d0 > d1)
                {
                    which = 0;
                    transStrand = 1;
                }
                else if (d0 < d1)
                {
                    which = 1;
                    transStrand = 2;
                }
                else
                {
                    transStrand = 3;
                    which = (qlen + d0) & 1;
                }

                if (which == 0)
                {
                    regs[i] = s0;
                    r2 = s20;
                }
                else
                {
                    regs[i] = s1;
                    r2 = s21;
                }

                if (regs[i].P is not null)
                    regs[i].P!.TransStrand = transStrand;
            }
            else
            {
                r2 = Align1(o, minDpMax, mi, qlen, qseq0, regs[i], nA, a, ez, o.SpliceForward, o.SpliceReverse);
                if (o.Splice && regs[i].P is not null)
                    regs[i].P!.TransStrand = o.SpliceForward ? 1 : 2;
            }

            if (r2 is not null && r2.Cnt > 0)
                regs.Insert(i + 1, r2);
            if (i > 0 && regs[i].SplitInv && !o.NoInversion)
            {
                var inv = Align1Inv(o, minDpMax, mi, qlen, qseq0, regs[i - 1], regs[i], ez);
                if (inv is not null)
                {
                    regs.Insert(i + 1, inv);
                    ++i;
                }
            }
        }

        FilterRegs(o, minDpMax, qlen, regs);
        if (!o.ShortRead && !o.SplitPrefix && qlen >= o.RankMinLength)
        {
            UpdateDpMax(qlen, regs, o.RankFraction, o.Scoring.Match, o.Scoring.Mismatch);
            FilterRegs(o, minDpMax, qlen, regs);
        }

        HitSort(regs, o.AltDrop);
        return regs;
    }
}
