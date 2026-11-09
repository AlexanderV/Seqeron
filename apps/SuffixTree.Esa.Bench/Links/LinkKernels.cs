#pragma warning disable S2325 // stubs of the red run: the bodies are written after the tests
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

/// <summary>The suffix-link kernels of prototype P1 (work package 1b.4a).</summary>
internal static class LinkKernels
{
    /// <summary>The kernel as a delegate; the delegate indirection is the only overhead (a few ns against links of 100 ns and more).</summary>
    internal static LinkFunc<TCnt> Get<TCnt>(LinkKernel kernel, EndpointSource source)
        where TCnt : struct, IStepCounter => throw new NotImplementedException();

    /// <summary>The link of <paramref name="node"/> by <paramref name="kernel"/> (no step limit).</summary>
    internal static LinkTarget Link(LinkKernel kernel, EndpointSource source, LinkContext ctx, in LinkNode node) => throw new NotImplementedException();

    /// <summary>The link with a step budget; false (and <paramref name="steps"/> above the budget) when the kernel gave up.</summary>
    internal static bool TryLink(LinkKernel kernel, EndpointSource source, LinkContext ctx, in LinkNode node, long budget, out LinkTarget target, out long steps)
        => throw new NotImplementedException();
}
