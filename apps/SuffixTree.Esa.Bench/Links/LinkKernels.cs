using System.Runtime.CompilerServices;

namespace SuffixTree.Esa.Bench.Links;

/// <summary>Which of the interval searches of ADR-001 section 5 finds the link target.</summary>
internal enum LinkKernel
{
    /// <summary>(a) LCP outward scan from <c>[Psi(l), Psi(r)]</c> while <c>LCP &gt;= depth - 1</c>.</summary>
    Scan,

    /// <summary>
    /// (b) child-table descent along the rank path of <c>Psi(l)</c> from the root until the node of value <c>depth - 1</c> (the layout stores no parent pointers, so a parent
    /// climb is not available; see the Log of PROGRESS-S1.md).
    /// </summary>
    Descent,

    /// <summary>(c) as in the ADR: <c>p = RMQ(a+1..b)</c>, <c>L = PSV_strict(p)</c>, <c>R + 1 = NSV_strict(p)</c> on the block-minima pyramid.</summary>
    RmqPsvNsv,

    /// <summary>(c2) the value-parameterised form of (c): <c>L = PrevLess(a, t)</c>, <c>R + 1 = NextLess(b + 1, t)</c>; no range-minimum query.</summary>
    PsvNsv,
}

/// <summary>Where the endpoints <c>Psi(l)</c>, <c>Psi(r)</c> come from.</summary>
internal enum EndpointSource
{
    /// <summary><c>ISA[SA[l] + 1]</c> and <c>ISA[SA[r] + 1]</c> (the node handle already carries <c>SA[l]</c>).</summary>
    Isa,

    /// <summary><c>Psi[l]</c>, <c>Psi[r]</c> (one independent read per endpoint).</summary>
    Psi,
}

/// <summary>An internal node (lcp-interval) as the navigator handle will carry it: rank interval, string depth and <c>Pos = SA[Lb]</c>.</summary>
internal readonly record struct LinkNode(int Lb, int Rb, int Depth, int Pos);

/// <summary>The node a suffix link leads to: its rank interval and the two facts a child-table walk needs (parent depth; whether it is the last child). The root is <c>(0, n, -1, true)</c>.</summary>
internal readonly record struct LinkTarget(int Lb, int Rb, int ParentDepth, bool IsLast);

/// <summary>A link kernel as a delegate; false when the step budget of the counter was exceeded (the target is then undefined).</summary>
internal delegate bool LinkFunc<TCnt>(LinkContext ctx, in LinkNode node, ref TCnt counter, out LinkTarget target)
    where TCnt : struct, IStepCounter;

/// <summary>The endpoint source as a type, so that a kernel compiled for it contains only its reads.</summary>
internal interface IEndpointSource
{
    /// <summary><c>Psi(Lb)</c>.</summary>
    static abstract int Left(LinkContext c, in LinkNode v);

    /// <summary><c>Psi(Rb)</c>.</summary>
    static abstract int Right(LinkContext c, in LinkNode v);

    /// <summary>Throws <see cref="InvalidOperationException"/> when the context lacks the array.</summary>
    static abstract void Require(LinkContext c);
}

/// <summary><c>ISA[SA[l] + 1]</c> and <c>ISA[SA[r] + 1]</c>: the left endpoint costs one read (the handle carries <c>SA[l]</c>), the right one two dependent reads.</summary>
internal readonly struct IsaSource : IEndpointSource
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int Left(LinkContext c, in LinkNode v) => c.Isa![v.Pos + 1];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int Right(LinkContext c, in LinkNode v) => c.Isa![c.Sa[v.Rb] + 1];

    public static void Require(LinkContext c)
    {
        if (c.Isa is null) throw new InvalidOperationException("this endpoint source needs ISA");
    }
}

/// <summary><c>Psi[l]</c> and <c>Psi[r]</c>: one independent read per endpoint.</summary>
internal readonly struct PsiSource : IEndpointSource
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int Left(LinkContext c, in LinkNode v) => c.Psi![v.Lb];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int Right(LinkContext c, in LinkNode v) => c.Psi![v.Rb];

    public static void Require(LinkContext c)
    {
        if (c.Psi is null) throw new InvalidOperationException("this endpoint source needs Psi");
    }
}

/// <summary>The suffix-link kernels of prototype P1 (work package 1b.4a).</summary>
internal static class LinkKernels
{
    private static class Funcs<TCnt>
        where TCnt : struct, IStepCounter
    {
        internal static readonly LinkFunc<TCnt>[,] Table = Build();

        private static LinkFunc<TCnt>[,] Build()
        {
            var t = new LinkFunc<TCnt>[4, 2];
            t[(int)LinkKernel.Scan, (int)EndpointSource.Isa] = Scan<IsaSource, TCnt>;
            t[(int)LinkKernel.Scan, (int)EndpointSource.Psi] = Scan<PsiSource, TCnt>;
            t[(int)LinkKernel.Descent, (int)EndpointSource.Isa] = Descent<IsaSource, TCnt>;
            t[(int)LinkKernel.Descent, (int)EndpointSource.Psi] = Descent<PsiSource, TCnt>;
            t[(int)LinkKernel.RmqPsvNsv, (int)EndpointSource.Isa] = RmqPsvNsv<IsaSource, TCnt>;
            t[(int)LinkKernel.RmqPsvNsv, (int)EndpointSource.Psi] = RmqPsvNsv<PsiSource, TCnt>;
            t[(int)LinkKernel.PsvNsv, (int)EndpointSource.Isa] = PsvNsv<IsaSource, TCnt>;
            t[(int)LinkKernel.PsvNsv, (int)EndpointSource.Psi] = PsvNsv<PsiSource, TCnt>;
            return t;
        }
    }

    /// <summary>The kernel as a delegate; the delegate indirection is the only overhead (a few ns against links of 100 ns and more).</summary>
    internal static LinkFunc<TCnt> Get<TCnt>(LinkKernel kernel, EndpointSource source)
        where TCnt : struct, IStepCounter
    {
        if (!Enum.IsDefined(kernel)) throw new ArgumentOutOfRangeException(nameof(kernel));
        if (!Enum.IsDefined(source)) throw new ArgumentOutOfRangeException(nameof(source));
        return Funcs<TCnt>.Table[(int)kernel, (int)source];
    }

    /// <summary>The link of <paramref name="node"/> by <paramref name="kernel"/> (no step limit).</summary>
    internal static LinkTarget Link(LinkKernel kernel, EndpointSource source, LinkContext ctx, in LinkNode node)
    {
        var none = new NoCount();
        Get<NoCount>(kernel, source)(ctx, node, ref none, out LinkTarget target);
        return target;
    }

    /// <summary>The link with a step budget; false (and <paramref name="steps"/> above the budget) when the kernel gave up.</summary>
    internal static bool TryLink(LinkKernel kernel, EndpointSource source, LinkContext ctx, in LinkNode node, long budget, out LinkTarget target, out long steps)
    {
        var counter = new BudgetSteps(budget);
        bool ok = Get<BudgetSteps>(kernel, source)(ctx, node, ref counter, out target);
        steps = counter.Steps;
        return ok;
    }

    private static LinkTarget Root(LinkContext ctx) => new(0, ctx.N, -1, true);

    private static void RequireHelper(LinkContext ctx)
    {
        if (ctx.Helper is null) throw new InvalidOperationException("this kernel needs the LCP range helper");
    }

    /// <summary>(a) outward scan.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static bool Scan<TSrc, TCnt>(LinkContext ctx, in LinkNode v, ref TCnt cnt, out LinkTarget target)
        where TSrc : struct, IEndpointSource
        where TCnt : struct, IStepCounter
    {
        TSrc.Require(ctx);
        int t = v.Depth - 1;
        if (t == 0)
        {
            target = Root(ctx);
            return true;
        }

        int a = TSrc.Left(ctx, v), b = TSrc.Right(ctx, v);
        CompactLcp lcp = ctx.Lcp;
        int n = ctx.N;
        int l = a;
        while (true)
        {
            cnt.Add(1);
            if (lcp[l] < t) break;
            l--;
            if (cnt.Exceeded)
            {
                target = default;
                return false;
            }
        }
        int r = b + 1;
        while (r <= n)
        {
            cnt.Add(1);
            if (lcp[r] < t) break;
            r++;
            if (cnt.Exceeded)
            {
                target = default;
                return false;
            }
        }
        target = Boundary(lcp, n, l, r);
        return !cnt.Exceeded;
    }

    /// <summary>The target <c>[l .. r - 1]</c> with the two LCP values around it (cached by the search that just read them).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static LinkTarget Boundary(CompactLcp lcp, int n, int l, int r)
    {
        int left = lcp[l];
        int right = r <= n ? lcp[r] : -1;
        return new LinkTarget(l, r - 1, Math.Max(left, right), left > right);
    }

    /// <summary>(b) child-table descent along the rank path of <c>Psi(l)</c>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static bool Descent<TSrc, TCnt>(LinkContext ctx, in LinkNode v, ref TCnt cnt, out LinkTarget target)
        where TSrc : struct, IEndpointSource
        where TCnt : struct, IStepCounter
    {
        TSrc.Require(ctx);
        int t = v.Depth - 1;
        if (t == 0)
        {
            target = Root(ctx);
            return true;
        }

        int a = TSrc.Left(ctx, v);
        TaggedChild child = ctx.Child;
        CompactLcp lcp = ctx.Lcp;
        int lb = 0, rb = ctx.N, ell = 0;
        bool last = true;
        while (true)
        {
            // the child of [lb..rb] that holds rank a: the last child that starts at or before a
            int boundary = child.FirstLIndex(lb, rb, last);
            cnt.Add(1);
            int cl = lb;
            while (boundary >= 0 && boundary <= a)
            {
                cl = boundary;
                boundary = child.NextLIndex(boundary, rb);
                cnt.Add(1);
            }
            int cr = boundary >= 0 ? boundary - 1 : rb;
            bool cLast = boundary < 0;
            if (cl == cr) throw new InvalidOperationException($"the descent to rank {a} reached a leaf above depth {t}: no node of that value on the path (identity F1 violated)");

            int cell = lcp[child.FirstLIndex(cl, cr, cLast)];
            cnt.Add(1);
            if (cell == t)
            {
                target = new LinkTarget(cl, cr, ell, cLast);
                return true;
            }
            if (cell > t) throw new InvalidOperationException($"the descent to rank {a} passed depth {t} (child value {cell}): identity F1 violated");
            if (cnt.Exceeded)
            {
                target = default;
                return false;
            }
            lb = cl;
            rb = cr;
            ell = cell;
            last = cLast;
        }
    }

    /// <summary>(c) range minimum with its position, then previous and next smaller value of that position.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static bool RmqPsvNsv<TSrc, TCnt>(LinkContext ctx, in LinkNode v, ref TCnt cnt, out LinkTarget target)
        where TSrc : struct, IEndpointSource
        where TCnt : struct, IStepCounter
    {
        TSrc.Require(ctx);
        RequireHelper(ctx);
        int t = v.Depth - 1;
        if (t == 0)
        {
            target = Root(ctx);
            return true;
        }

        LcpRangeHelper h = ctx.Helper!;
        int a = TSrc.Left(ctx, v), b = TSrc.Right(ctx, v);
        int p = h.ArgMin(a + 1, b, ref cnt);
        cnt.Add(1);
        int value = ctx.Lcp[p];
        int l = h.PrevLess(p - 1, value, ref cnt);
        int r = h.NextLess(p + 1, value, ref cnt);
        target = Boundary(ctx.Lcp, ctx.N, l, r);
        return !cnt.Exceeded;
    }

    /// <summary>(c2) previous / next smaller than the known value <c>t</c> from the two endpoints; no range-minimum query.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static bool PsvNsv<TSrc, TCnt>(LinkContext ctx, in LinkNode v, ref TCnt cnt, out LinkTarget target)
        where TSrc : struct, IEndpointSource
        where TCnt : struct, IStepCounter
    {
        TSrc.Require(ctx);
        RequireHelper(ctx);
        int t = v.Depth - 1;
        if (t == 0)
        {
            target = Root(ctx);
            return true;
        }

        LcpRangeHelper h = ctx.Helper!;
        int a = TSrc.Left(ctx, v), b = TSrc.Right(ctx, v);
        int l = h.PrevLess(a, t, ref cnt);
        int r = h.NextLess(b + 1, t, ref cnt);
        target = Boundary(ctx.Lcp, ctx.N, l, r);
        return !cnt.Exceeded;
    }
}
