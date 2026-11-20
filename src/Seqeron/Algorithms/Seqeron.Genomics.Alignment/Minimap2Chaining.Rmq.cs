namespace Seqeron.Genomics.Alignment;

public static partial class Minimap2Chaining
{
    // ------------------------------------------------------------------------------------------------ krmq.h: an AVL tree with range-minimum queries (the structure of mg_lchain_rmq)
    // Ported exactly (insert, erase, rotations, rmq, interval, iterator): the node that a range-minimum query returns among equal priorities depends on the shape of the tree, which depends on the history of insertions and
    // erasures, so the tree itself is the program's.

    private sealed class LcElem
    {
        public int Y;
        public long I;
        public double Pri;
        public readonly LcElem?[] P = new LcElem?[2];
        public LcElem S = null!;
        public sbyte Balance;
        public uint Size;
    }

    private const int RmqMaxDepth = 64;

    private static int LcCmp(LcElem a, LcElem b)
    {
        if (a.Y != b.Y)
            return a.Y < b.Y ? -1 : 1;
        return a.I.CompareTo(b.I);
    }

    private static bool LcLt2(LcElem a, LcElem b) => a.Pri < b.Pri;

    private static uint TreeSize(LcElem? p) => p?.Size ?? 0;

    private static uint SizeChild(LcElem q, int i) => q.P[i]?.Size ?? 0;

    private static void UpdateMin(LcElem p, LcElem? q, LcElem? r)
    {
        p.S = q is null || LcLt2(p, q.S) ? p : q.S;
        p.S = r is null || LcLt2(p.S, r.S) ? p.S : r.S;
    }

    // one rotation: (a,(b,c)q)p => ((a,b)p,c)q
    private static LcElem Rotate1(LcElem p, int dir)
    {
        int opp = 1 - dir;
        var q = p.P[opp]!;
        var s = p.S;
        uint sizeP = p.Size;
        p.Size -= q.Size - SizeChild(q, dir);
        q.Size = sizeP;
        UpdateMin(p, p.P[dir], q.P[dir]);
        q.S = s;
        p.P[opp] = q.P[dir];
        q.P[dir] = p;
        return q;
    }

    // two consecutive rotations: (a,((b,c)r,d)q)p => ((a,b)p,(c,d)q)r
    private static LcElem Rotate2(LcElem p, int dir)
    {
        int opp = 1 - dir;
        var q = p.P[opp]!;
        var r = q.P[dir]!;
        var s = p.S;
        uint sizeXDir = SizeChild(r, dir);
        r.Size = p.Size;
        p.Size -= q.Size - sizeXDir;
        q.Size -= sizeXDir + 1;
        UpdateMin(p, p.P[dir], r.P[dir]);
        UpdateMin(q, q.P[opp], r.P[opp]);
        r.S = s;
        p.P[opp] = r.P[dir];
        r.P[dir] = p;
        q.P[dir] = r.P[opp];
        r.P[opp] = q;
        int b1 = dir == 0 ? +1 : -1;
        if (r.Balance == b1)
        {
            q.Balance = 0;
            p.Balance = (sbyte)-b1;
        }
        else if (r.Balance == 0)
        {
            q.Balance = p.Balance = 0;
        }
        else
        {
            q.Balance = (sbyte)b1;
            p.Balance = 0;
        }

        r.Balance = 0;
        return r;
    }

    private static void Insert(ref LcElem? root, LcElem x)
    {
        var stack = new byte[RmqMaxDepth];
        var path = new LcElem[RmqMaxDepth];
        LcElem? bp = root, bq = null, p, q;
        int which = 0, top = 0, pathLen = 0;
        for (p = bp, q = bq; p is not null; q = p, p = p.P[which])
        {
            int cmp = LcCmp(x, p);
            if (cmp == 0)
                return;
            if (p.Balance != 0)
            {
                bq = q;
                bp = p;
                top = 0;
            }

            which = cmp > 0 ? 1 : 0;
            stack[top++] = (byte)which;
            path[pathLen++] = p;
        }

        x.Balance = 0;
        x.Size = 1;
        x.P[0] = x.P[1] = null;
        x.S = x;
        if (q is null)
            root = x;
        else
            q.P[which] = x;
        if (bp is null)
            return;
        for (int i = 0; i < pathLen; ++i)
            ++path[i].Size;
        for (int i = pathLen - 1; i >= 0; --i)
        {
            UpdateMin(path[i], path[i].P[0], path[i].P[1]);
            if (path[i].S != x)
                break;
        }

        top = 0;
        for (p = bp; p != x; p = p!.P[stack[top]], ++top)
        {
            if (stack[top] == 0)
                --p!.Balance;
            else
                ++p!.Balance;
        }

        if (bp.Balance is > -2 and < 2)
            return;
        which = bp.Balance < 0 ? 1 : 0;
        int b1 = which == 0 ? +1 : -1;
        q = bp.P[1 - which]!;
        LcElem r;
        if (q.Balance == b1)
        {
            r = Rotate1(bp, which);
            q.Balance = bp.Balance = 0;
        }
        else
        {
            r = Rotate2(bp, which);
        }

        if (bq is null)
            root = r;
        else
            bq.P[bp != bq.P[0] ? 1 : 0] = r;
        return;
    }

    private static LcElem? Find(LcElem? root, LcElem x)
    {
        var p = root;
        while (p is not null)
        {
            int cmp = LcCmp(x, p);
            if (cmp < 0)
                p = p.P[0];
            else if (cmp > 0)
                p = p.P[1];
            else
                break;
        }

        return p;
    }

    private static void Interval(LcElem? root, LcElem x, out LcElem? lower, out LcElem? upper)
    {
        var p = root;
        LcElem? l = null, u = null;
        while (p is not null)
        {
            int cmp = LcCmp(x, p);
            if (cmp < 0)
            {
                u = p;
                p = p.P[0];
            }
            else if (cmp > 0)
            {
                l = p;
                p = p.P[1];
            }
            else
            {
                l = u = p;
                break;
            }
        }

        lower = l;
        upper = u;
    }

    private static LcElem? Rmq(LcElem? root, LcElem lo, LcElem up)
    {
        if (root is null)
            return null;
        var path = new LcElem[2][] { new LcElem[RmqMaxDepth], new LcElem[RmqMaxDepth] };
        var pcmp = new int[2][] { new int[RmqMaxDepth], new int[RmqMaxDepth] };
        var plen = new int[2];
        var p = root;
        while (p is not null)
        {
            int cmp = LcCmp(lo, p);
            path[0][plen[0]] = p;
            pcmp[0][plen[0]++] = cmp;
            if (cmp < 0)
                p = p.P[0];
            else if (cmp > 0)
                p = p.P[1];
            else
                break;
        }

        p = root;
        while (p is not null)
        {
            int cmp = LcCmp(up, p);
            path[1][plen[1]] = p;
            pcmp[1][plen[1]++] = cmp;
            if (cmp < 0)
                p = p.P[0];
            else if (cmp > 0)
                p = p.P[1];
            else
                break;
        }

        int i;
        for (i = 0; i < plen[0] && i < plen[1]; ++i) // find the LCA
            if (path[0][i] == path[1][i] && pcmp[0][i] <= 0 && pcmp[1][i] >= 0)
                break;
        if (i == plen[0] || i == plen[1])
            return null; // no elements in the closed interval
        int lca = i;
        var min = path[0][lca];
        for (i = lca + 1; i < plen[0]; ++i)
        {
            if (pcmp[0][i] <= 0)
            {
                if (LcLt2(path[0][i], min))
                    min = path[0][i];
                var right = path[0][i].P[1];
                if (right is not null && LcLt2(right.S, min))
                    min = right.S;
            }
        }

        for (i = lca + 1; i < plen[1]; ++i)
        {
            if (pcmp[1][i] >= 0)
            {
                if (LcLt2(path[1][i], min))
                    min = path[1][i];
                var left = path[1][i].P[0];
                if (left is not null && LcLt2(left.S, min))
                    min = left.S;
            }
        }

        return min;
    }

    private static void Erase(ref LcElem? root, LcElem x)
    {
        var path = new LcElem[RmqMaxDepth];
        var dir = new byte[RmqMaxDepth];
        int d = 0;
        var fake = new LcElem();
        fake.P[0] = root;
        fake.P[1] = null;
        LcElem? p = fake;
        int cmp = -1;
        while (cmp != 0)
        {
            int which = cmp > 0 ? 1 : 0;
            dir[d] = (byte)which;
            path[d++] = p!;
            p = p!.P[which];
            if (p is null)
                return;
            cmp = LcCmp(x, p);
        }

        for (int i = 1; i < d; ++i)
            --path[i].Size;
        if (p!.P[1] is null)
        {
            // ((1,.)2,3)4 => (1,3)4; p=2
            path[d - 1].P[dir[d - 1]] = p.P[0];
        }
        else
        {
            var q = p.P[1]!;
            if (q.P[0] is null)
            {
                // ((1,2)3,4)5 => ((1)2,4)5; p=3,q=2
                q.P[0] = p.P[0];
                q.Balance = p.Balance;
                path[d - 1].P[dir[d - 1]] = q;
                path[d] = q;
                dir[d++] = 1;
                q.Size = p.Size - 1;
            }
            else
            {
                // ((1,((.,2)3,4)5)6,7)8 => ((1,(2,4)5)3,7)8; p=6
                LcElem r;
                int e = d++;
                while (true)
                {
                    dir[d] = 0;
                    path[d++] = q;
                    r = q.P[0]!;
                    if (r.P[0] is null)
                        break;
                    q = r;
                }

                r.P[0] = p.P[0];
                q.P[0] = r.P[1];
                r.P[1] = p.P[1];
                r.Balance = p.Balance;
                path[e - 1].P[dir[e - 1]] = r;
                path[e] = r;
                dir[e] = 1;
                for (int i = e + 1; i < d; ++i)
                    --path[i].Size;
                r.Size = p.Size - 1;
            }
        }

        for (int i = d - 1; i >= 0; --i)
            UpdateMin(path[i], path[i].P[0], path[i].P[1]);
        while (--d > 0)
        {
            var q = path[d];
            int b1 = 1, b2 = 2;
            int which = dir[d], other = 1 - which;
            if (which != 0)
            {
                b1 = -b1;
                b2 = -b2;
            }

            q.Balance += (sbyte)b1;
            if (q.Balance == b1)
                break;
            if (q.Balance == b2)
            {
                var r = q.P[other]!;
                if (r.Balance == -b1)
                {
                    path[d - 1].P[dir[d - 1]] = Rotate2(q, which);
                }
                else
                {
                    path[d - 1].P[dir[d - 1]] = Rotate1(q, which);
                    if (r.Balance == 0)
                    {
                        r.Balance = (sbyte)-b1;
                        q.Balance = (sbyte)b1;
                        break;
                    }

                    r.Balance = q.Balance = 0;
                }
            }
        }

        root = fake.P[0];
    }

    private struct RmqIterator
    {
        public LcElem[] Stack;
        public int Top; // -1 = before the first element

        public readonly LcElem? At => Top < 0 ? null : Stack[Top];
    }

    private static void IteratorFind(LcElem? root, LcElem x, out RmqIterator itr)
    {
        itr = new RmqIterator { Stack = new LcElem[RmqMaxDepth], Top = -1 };
        var p = root;
        while (p is not null)
        {
            itr.Stack[++itr.Top] = p;
            int cmp = LcCmp(x, p);
            if (cmp < 0)
                p = p.P[0];
            else if (cmp > 0)
                p = p.P[1];
            else
                break;
        }
    }

    // krmq_itr_next_bidir with dir = 0: the previous element
    private static bool IteratorPrev(ref RmqIterator itr)
    {
        if (itr.Top < 0)
            return false;
        const int dir = 0;
        var p = itr.Stack[itr.Top].P[dir];
        if (p is not null)
        {
            for (; p is not null; p = p.P[1 - dir])
                itr.Stack[++itr.Top] = p;
            return true;
        }

        do
        {
            p = itr.Stack[itr.Top--];
        }
        while (itr.Top >= 0 && p == itr.Stack[itr.Top].P[dir]);
        return itr.Top >= 0;
    }

    // comput_sc_simple
    private static int ScoreSimple(Mm128 ai, Mm128 aj, float chnPenGap, float chnPenSkip, out bool exact, out int width)
    {
        int dq = (int)ai.Y - (int)aj.Y;
        int dr = (int)(ai.X - aj.X);
        width = dr > dq ? dr - dq : dq - dr;
        int dd = width;
        int dg = dr < dq ? dr : dq;
        int qSpan = (int)(aj.Y >> 32 & 0xff);
        int sc = qSpan < dg ? qSpan : dg;
        exact = dd == 0 && dg <= qSpan;
        if (dd != 0 || dq > qSpan)
        {
            float linPen = chnPenGap * dd + chnPenSkip * dg;
            float logPen = dd >= 1 ? Log2(dd + 1) : 0.0f;
            sc = unchecked(sc - CInt(linPen + .5f * logPen));
        }

        return sc;
    }

    /// <summary><c>mg_lchain_rmq</c> (lchain.c): the chaining with a range-minimum-query tree over the anchors in range; the chains as <c>compact_a</c> leaves them. <paramref name="a"/> must be sorted by <c>x</c>.</summary>
    private static IReadOnlyList<MinimapChain> RmqChain(Mm128[] a, int maxDist, int maxDistInner, int bw, int maxChnSkip, int capRmqSize, int minCnt, int minSc, float chnPenGap, float chnPenSkip)
    {
        int n = a.Length;
        if (n == 0)
            return [];
        int maxDrop = bw;
        if (maxDist < bw)
            maxDist = bw;
        if (maxDistInner <= 0 || maxDistInner >= maxDist)
            maxDistInner = 0;
        var f = new int[n];
        var t = new int[n];
        var p = new long[n];
        LcElem? root = null, rootInner = null;
        long st = 0, stInner = 0, i0 = 0;
        for (long i = 0; i < n; ++i)
        {
            if ((i & 63) == 0)
                AlignmentCancellation.ThrowIfRequested();
            long maxJ = -1;
            int maxF = (int)(a[i].Y >> 32 & 0xff);
            // add in-range anchors
            if (i0 < i && a[i0].X != a[i].X)
            {
                for (long j = i0; j < i; ++j)
                {
                    var q = new LcElem { Y = (int)a[j].Y, I = j, Pri = -(f[j] + 0.5 * chnPenGap * ((int)a[j].X + (int)a[j].Y)) };
                    Insert(ref root, q);
                    if (maxDistInner > 0)
                        Insert(ref rootInner, new LcElem { Y = q.Y, I = q.I, Pri = q.Pri });
                }

                i0 = i;
            }

            // get rid of active chains out of range
            while (st < i && (a[i].X >> 32 != a[st].X >> 32 || a[i].X > a[st].X + (ulong)maxDist || TreeSize(root) > (uint)capRmqSize))
            {
                var s = new LcElem { Y = (int)a[st].Y, I = st };
                var found = Find(root, s);
                if (found is not null)
                    Erase(ref root, found);
                ++st;
            }

            if (maxDistInner > 0)
            {
                while (stInner < i && (a[i].X >> 32 != a[stInner].X >> 32 || a[i].X > a[stInner].X + (ulong)maxDistInner || TreeSize(rootInner) > (uint)capRmqSize))
                {
                    var s = new LcElem { Y = (int)a[stInner].Y, I = stInner };
                    var found = Find(rootInner, s);
                    if (found is not null)
                        Erase(ref rootInner, found);
                    ++stInner;
                }
            }

            // RMQ
            var lo = new LcElem { I = int.MaxValue, Y = (int)a[i].Y - maxDist };
            var hi = new LcElem { I = 0, Y = (int)a[i].Y };
            var qr = Rmq(root, lo, hi);
            if (qr is not null)
            {
                int nSkip = 0;
                long j = qr.I;
                int sc = f[j] + ScoreSimple(a[i], a[j], chnPenGap, chnPenSkip, out bool exact, out int width);
                if (width <= bw && sc > maxF)
                {
                    maxF = sc;
                    maxJ = j;
                }

                if (!exact && rootInner is not null && (int)a[i].Y > 0)
                {
                    var s = new LcElem { Y = (int)a[i].Y - 1, I = n };
                    Interval(rootInner, s, out var lower, out _);
                    if (lower is not null)
                    {
                        IteratorFind(rootInner, lower, out var itr);
                        LcElem? qq;
                        while ((qq = itr.At) is not null)
                        {
                            if (qq.Y < (int)a[i].Y - maxDistInner)
                                break;
                            j = qq.I;
                            sc = f[j] + ScoreSimple(a[i], a[j], chnPenGap, chnPenSkip, out _, out width);
                            if (width <= bw)
                            {
                                if (sc > maxF)
                                {
                                    maxF = sc;
                                    maxJ = j;
                                    if (nSkip > 0)
                                        --nSkip;
                                }
                                else if (t[j] == (int)i && ++nSkip > maxChnSkip)
                                {
                                    break;
                                }

                                if (p[j] >= 0)
                                    t[p[j]] = (int)i;
                            }

                            if (!IteratorPrev(ref itr))
                                break;
                        }
                    }
                }
            }

            f[i] = maxF;
            p[i] = maxJ;
        }

        return Backtrack(n, f, p, t, minCnt, minSc, maxDrop, a.Select(m => m.X).ToArray());
    }
}
