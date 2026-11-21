// The control flow of the ported C functions is kept literally (loop variables moved inside the body, empty search loops, the C compound conditions).
#pragma warning disable S127, S108, S1994, S3358, S1066, S1854

namespace Seqeron.Genomics.Alignment;

// pe.c and the multi-segment part of hit.c of minimap2 2.26: the primary / secondary selection of a fragment chain (mm_select_sub_multi), the split of a fragment chain into the chains of its segments (mm_seg_gen), the pairing of the
// hits of two reads under -c (mm_pair) and the through-read flag (mm_set_pe_thru).
public static partial class Minimap2Chaining
{
    private const ulong SeedSegMask = 0xffUL << 48;
    private const int SeedSegShift = 48;

    // mm_select_sub_multi
    private static void SelectSubMulti(float priRatio, float pri1, float pri2, int maxGapRef, int minDiff, int bestN, int nSegs, int[] qlens, List<Reg> r)
    {
        if (!(priRatio > 0.0f) || r.Count <= 0)
            return;
        int n = r.Count, k = 0, n2nd = 0;
        int maxDist = nSegs == 2 ? qlens[0] + qlens[1] + maxGapRef : 0;
        for (int i = 0; i < n; ++i)
        {
            bool toKeep = false;
            if (r[i].Parent == i)
            {
                toKeep = true;
            }
            else if (r[i].Score + minDiff >= r[r[i].Parent].Score)
            {
                toKeep = true;
            }
            else
            {
                var p = r[r[i].Parent];
                var q = r[i];
                if (p.Rev == q.Rev && p.Rid == q.Rid && q.Re - p.Rs < maxDist && p.Re - q.Rs < maxDist)
                {
                    if (q.Score >= p.Score * pri1)
                        toKeep = true;
                }
                else
                {
                    bool isParBoth = nSegs == 2 && p.Qs < qlens[0] && p.Qe > qlens[0];
                    bool isChiBoth = nSegs == 2 && q.Qs < qlens[0] && q.Qe > qlens[0];
                    if (isChiBoth || isChiBoth == isParBoth)
                    {
                        if (q.Score >= p.Score * priRatio)
                            toKeep = true;
                    }
                    else if (q.Score >= p.Score * pri2)
                    {
                        toKeep = true;
                    }
                }
            }

            if (toKeep && r[i].Parent != i && n2nd++ >= bestN)
                toKeep = false;
            if (toKeep)
                r[k++] = r[i];
        }

        if (k != n)
        {
            r.RemoveRange(k, n - k);
            SyncRegs(r);
        }
    }

    // mm_set_pe_thru
    private static void SetPeThru(int[] qlens, List<Reg>[] regs)
    {
        var nPri = new int[2];
        var pri = new[] { -1, -1 };
        for (int s = 0; s < 2; ++s)
            for (int i = 0; i < regs[s].Count; ++i)
                if (regs[s][i].Id == regs[s][i].Parent)
                {
                    ++nPri[s];
                    pri[s] = i;
                }

        if (nPri[0] == 1 && nPri[1] == 1)
        {
            var p = regs[0][pri[0]];
            var q = regs[1][pri[1]];
            if (p.Rid == q.Rid && p.Rev == q.Rev && Math.Abs(p.Rs - q.Rs) < 3 && Math.Abs(p.Re - q.Re) < 3
                && ((p.Qs == 0 && qlens[1] - q.Qe == 0) || (q.Qs == 0 && qlens[0] - p.Qe == 0)))
            {
                p.PeThru = q.PeThru = true;
            }
        }
    }

    // mm_pair: pair the hits of the two reads (runs under -c only, it needs the DP maxima); sets proper_frag, lifts the paired hits to primary and raises the mapping qualities
    private static void Pair(int maxGapRef, int peBonus, int subDiff, int matchSc, int[] qlens, List<Reg>[] regs)
    {
        var items = new List<(int S, bool Rev, Reg R)>();
        int dpThres = 0, segs = 0;
        var keys = new List<Mm128>();
        for (int s = 0; s < 2; ++s)
        {
            int max = 0;
            for (int i = 0; i < regs[s].Count; ++i)
            {
                var r = regs[s][i];
                items.Add((s, r.Rev, r));
                keys.Add(new Mm128(((ulong)(uint)r.Rid << 32) | ((ulong)(uint)r.Rs << 1) | (r.Rev ? (uint)(s ^ 1) : (uint)s), (ulong)(items.Count - 1)));
                max = max > r.P!.DpMax ? max : r.P.DpMax;
                segs |= 1 << s;
            }

            dpThres += max;
        }

        if (segs != 3)
            return;
        dpThres -= peBonus;
        if (dpThres < 0)
            dpThres = 0;
        var karr = keys.ToArray();
        RadixSort(karr);
        var a = new List<(int S, bool Rev, Reg R, ulong Key)>(karr.Length);
        foreach (var kv in karr)
        {
            var it = items[(int)kv.Y];
            a.Add((it.S, it.Rev, it.R, kv.X));
        }

        int n = a.Count;
        long max2 = -1;
        var maxIdx = new[] { -1, -1 };
        var last = new[] { -1, -1 };
        var sc = new List<ulong>();
        for (int i = 0; i < n; ++i)
        {
            int revI = a[i].Rev ? 1 : 0;
            if ((a[i].Key & 1) != 0)
            {
                if (last[revI] < 0)
                    continue;
                var r = a[i].R;
                var q = a[last[revI]].R;
                if (r.Rid != q.Rid || r.Rs - q.Re > maxGapRef)
                    continue;
                for (int j = last[revI]; j >= 0; --j)
                {
                    if (a[j].Rev != a[i].Rev || a[j].S == a[i].S)
                        continue;
                    q = a[j].R;
                    if (r.Rid != q.Rid || r.Rs - q.Re > maxGapRef)
                        break;
                    if (r.P!.DpMax + q.P!.DpMax < dpThres)
                        continue;
                    long score = ((long)(r.P.DpMax + q.P.DpMax) << 32) | unchecked(r.Hash + q.Hash);
                    if (score > max2)
                    {
                        max2 = score;
                        maxIdx[a[j].S] = j;
                        maxIdx[a[i].S] = i;
                    }

                    sc.Add((ulong)score);
                }
            }
            else
            {
                last[revI] = i;
            }
        }

        if (sc.Count > 1)
            sc.Sort();
        if (sc.Count > 0 && max2 > 0)
        {
            int nSub = 0, mapqPe;
            var rr = new[] { a[maxIdx[0]].R, a[maxIdx[1]].R };
            rr[0].ProperFrag = rr[1].ProperFrag = true;
            for (int s = 0; s < 2; ++s)
            {
                if (rr[s].Id != rr[s].Parent)
                {
                    var p = regs[s][rr[s].Parent];
                    for (int i = 0; i < regs[s].Count; ++i)
                        if (regs[s][i].Parent == p.Id)
                            regs[s][i].Parent = rr[s].Id;
                    p.Mapq = 0;
                }

                if (!rr[s].SamPri)
                {
                    for (int i = 0; i < regs[s].Count; ++i)
                        regs[s][i].SamPri = false;
                    rr[s].SamPri = true;
                }
            }

            mapqPe = rr[0].Mapq > rr[1].Mapq ? rr[0].Mapq : rr[1].Mapq;
            ulong maxHi = (ulong)max2 >> 32;
            foreach (ulong v in sc)
                if ((v >> 32) + (ulong)subDiff >= maxHi)
                    ++nSub;
            if (sc.Count > 1)
            {
                int mapqPeAlt = CInt((6.02f * (maxHi - (sc[^2] >> 32)) / matchSc) - (4.343f * MathF.Log(nSub)));
                mapqPe = mapqPe < mapqPeAlt ? mapqPe : mapqPeAlt;
            }

            if (rr[0].Mapq < mapqPe)
                rr[0].Mapq = CInt((.2f * rr[0].Mapq) + (.8f * mapqPe) + .499f);
            if (rr[1].Mapq < mapqPe)
                rr[1].Mapq = CInt((.2f * rr[1].Mapq) + (.8f * mapqPe) + .499f);
            if (sc.Count == 1)
            {
                if (rr[0].Mapq < 2)
                    rr[0].Mapq = 2;
                if (rr[1].Mapq < 2)
                    rr[1].Mapq = 2;
            }
            else if (maxHi > sc[^2] >> 32)
            {
                if (rr[0].Mapq < 1)
                    rr[0].Mapq = 1;
                if (rr[1].Mapq < 1)
                    rr[1].Mapq = 1;
            }
        }

        SetPeThru(qlens, regs);
    }

    // mm_seg_gen: split the chains of a fragment into the chains of its segments; returns the hits of every segment and the anchor array of every segment
    private static (List<Reg>[] Regs, Mm128[][] A) SegGen(uint hash, int nSegs, int[] qlens, List<Reg> regs0, Mm128[] a)
    {
        var acc = new int[nSegs + 1];
        for (int s = 1; s < nSegs; ++s)
            acc[s] = acc[s - 1] + qlens[s - 1];
        int qlenSum = acc[nSegs - 1] + qlens[nSegs - 1];
        int n0 = regs0.Count;
        var u = new ulong[nSegs][];
        var nA = new int[nSegs];
        for (int s = 0; s < nSegs; ++s)
        {
            u[s] = new ulong[n0];
            for (int i = 0; i < n0; ++i)
                u[s][i] = (ulong)(uint)regs0[i].Score << 32;
        }

        for (int i = 0; i < n0; ++i)
        {
            var r = regs0[i];
            for (int j = 0; j < r.Cnt; ++j)
            {
                int sid = (int)((a[r.As + j].Y & SeedSegMask) >> SeedSegShift);
                ++u[sid][i];
                ++nA[sid];
            }
        }

        var nU = new int[nSegs];
        var sa = new Mm128[nSegs][];
        for (int s = 0; s < nSegs; ++s)
        {
            for (int i = 0; i < n0; ++i)
                if ((int)u[s][i] != 0)
                    u[s][nU[s]++] = u[s][i];
            sa[s] = new Mm128[nA[s]];
            nA[s] = 0;
        }

        for (int i = 0; i < n0; ++i)
        {
            var r = regs0[i];
            for (int j = 0; j < r.Cnt; ++j)
            {
                int sid = (int)((a[r.As + j].Y & SeedSegMask) >> SeedSegShift);
                var a1 = a[r.As + j];
                a1.Y -= (ulong)(long)((a1.X >> 63) != 0 ? qlenSum - (qlens[sid] + acc[sid]) : acc[sid]);
                sa[sid][nA[sid]++] = a1;
            }
        }

        var regs = new List<Reg>[nSegs];
        for (int s = 0; s < nSegs; ++s)
        {
            regs[s] = GenRegs(hash, qlens[s], [.. u[s].Take(nU[s])], sa[s], false);
            foreach (var r in regs[s])
            {
                r.SegSplit = true;
            }
        }

        return (regs, sa);
    }
}
