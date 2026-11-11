using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Contract;

/// <summary>How much of the tree <see cref="EsaContractChecker"/> visits: every point of every edge and every rank pair (thorough), or the ends of long edges and a seeded sample of pairs (sampled).</summary>
/// <param name="EdgePointCap">Edges with at most this many interior positions are checked at every position.</param>
/// <param name="PairSamples">Rank pairs for <c>Lcp</c> / <c>Lca</c>: all pairs when (n + 1)^2 does not exceed it, else this many seeded pairs.</param>
/// <param name="Seed">Seed of the sampling.</param>
internal sealed record CheckOptions(int EdgePointCap, int PairSamples, ulong Seed)
{
    /// <summary>Every locus, every pair.</summary>
    internal static CheckOptions Thorough { get; } = new(int.MaxValue, int.MaxValue, 1);

    /// <summary>Long edges at their ends and a few seeded interior positions, a seeded sample of pairs.</summary>
    internal static CheckOptions Sampled { get; } = new(10, 600, 1);
}

/// <summary>
/// Checks every operation of an <see cref="EsaSuffixIndex"/> (or of the facade, work package 1b.6: the checks use the contract interfaces only) on every kind of locus (the root, internal nodes, leaves including the terminator leaves and the leaf of the empty suffix,
/// positions inside an edge) against <see cref="IndexContractOracle"/>, which is derived from the definitions. The tree is walked from the oracle's side: for every node string the engine
/// is asked for the locus (<c>TryChild</c> / <c>TryExtendRight</c>) and the answer is compared with the oracle's locus, then every operation is applied to it. Nothing here reads the engine's arrays.
/// </summary>
internal static class EsaContractChecker
{
    private static readonly IndexLocus Sentinel = new(7, 7, 7, 7);

    private static void Eq<T>(T actual, T expected, Func<string> what)
    {
        if (!EqualityComparer<T>.Default.Equals(actual, expected)) Assert.Fail($"{what()}: expected {expected}, got {actual}");
    }

    private static void IsTrue(bool condition, Func<string> what)
    {
        if (!condition) Assert.Fail(what());
    }

    /// <summary>Runs every check on <paramref name="idx"/>, which must index <paramref name="o"/>'s text.</summary>
    internal static void CheckAll<TIdx>(TIdx idx, IndexContractOracle o, string context, CheckOptions options)
        where TIdx : ISuffixTreeNavigation, ISuffixArrayAccess
    {
        CheckMetadata(idx, o, context);
        CheckArrays(idx, o, context, options);
        CheckTree(idx, o, context, options);
    }

    private static int[] Alphabet(IndexContractOracle o)
    {
        var symbols = new List<int>(o.Alphabet());
        int added = 0;
        foreach (int candidate in new[] { 'z', 0xFFFF, 0, 'A', 0x100 })
            if (!symbols.Contains(candidate) && added++ < 2) symbols.Add(candidate);
        return [.. symbols];
    }

    private static string P(int[] w) => IndexContractOracle.Render(w);

    // ------------------------------------------------------------------------------------------------------------------------------

    private static void CheckMetadata<TIdx>(TIdx idx, IndexContractOracle o, string context)
        where TIdx : ISuffixTreeNavigation, ISuffixArrayAccess
    {
        Eq(idx.Engine, "esa", () => context + " Engine");
        Eq(idx.TextLength, (long)o.N, () => context + " TextLength");
        Eq(idx.Text.Length, (long)o.N, () => context + " Text.Length");
        for (int i = 0; i < o.N; i++) Eq(idx.Text.SymbolAt(i), (int)o.Text[i], () => $"{context} SymbolAt({i})");
        Eq(idx.Text.SymbolAt(o.N), -1, () => context + " SymbolAt(n)");
        Eq(idx.Root, o.LocusOf([]), () => context + " Root");
        Eq(idx.Root, new IndexLocus(0, o.N, 0, 0), () => context + " Root representation");
        Eq(idx.Count(idx.Root), (long)o.N, () => context + " Count(Root)");
    }

    private static void CheckArrays<TIdx>(TIdx idx, IndexContractOracle o, string context, CheckOptions options)
        where TIdx : ISuffixTreeNavigation, ISuffixArrayAccess
    {
        int n = o.N;
        for (int r = 0; r <= n; r++)
        {
            Eq(idx.SuffixAt(r), (long)o.Sa[r], () => $"{context} SuffixAt({r})");
            Eq(idx.LcpAt(r), (long)o.Lcp[r], () => $"{context} LcpAt({r})");
        }
        for (int p = 0; p <= n; p++) Eq(idx.RankOf(p), (long)o.Isa[p], () => $"{context} RankOf({p})");

        foreach (var (a, b) in Pairs(n, options))
        {
            Eq(idx.Lcp(a, b), (long)o.LcpOfRanks(a, b), () => $"{context} Lcp({a}, {b})");
            Eq(idx.Lcp(b, a), (long)o.LcpOfRanks(a, b), () => $"{context} Lcp({b}, {a})");
            IndexLocus expected = o.LocusOf(o.LcaString(a, b));
            Eq(idx.Lca(a, b), expected, () => $"{context} Lca({a}, {b})");
            Eq(idx.Lca(b, a), expected, () => $"{context} Lca({b}, {a})");
        }
    }

    private static IEnumerable<(int A, int B)> Pairs(int n, CheckOptions options)
    {
        long all = (long)(n + 1) * (n + 1);
        if (all <= options.PairSamples)
        {
            for (int a = 0; a <= n; a++)
                for (int b = 0; b <= n; b++) yield return (a, b);
            yield break;
        }
        var rng = new SplitMix64(options.Seed ^ (ulong)n);
        yield return (0, 0);
        yield return (0, n);
        yield return (n, n);
        yield return (1, n);
        for (int i = 0; i < options.PairSamples; i++) yield return (rng.NextInt(n + 1), rng.NextInt(n + 1));
        for (int i = 0; i < 40; i++)
        {
            int a = rng.NextInt(n + 1);
            yield return (a, Math.Min(n, a + 1));
            yield return (a, a);
        }
    }

    // ------------------------------------------------------------------------------------------------------------------------------

    private static void CheckTree<TIdx>(TIdx idx, IndexContractOracle o, string context, CheckOptions options)
        where TIdx : ISuffixTreeNavigation, ISuffixArrayAccess
    {
        int[] alphabet = Alphabet(o);
        var rng = new SplitMix64(options.Seed * 31 + (ulong)o.N);
        var stack = new Stack<(int[] W, IndexLocus L)>();
        stack.Push(([], idx.Root));
        while (stack.Count > 0)
        {
            var (w, locus) = stack.Pop();
            CheckNode(idx, o, w, locus, alphabet, context);
            foreach (var (key, childString) in o.Children(w))
            {
                bool found = idx.TryChild(locus, key, out IndexLocus child);
                IsTrue(found, () => $"{context} TryChild({P(w)}, {key}) must find {P(childString)}");
                Eq(child, o.LocusOf(childString), () => $"{context} TryChild({P(w)}, {key})");
                CheckEdge(idx, o, w, locus, childString, child, alphabet, context, options, rng);
                stack.Push((childString, child));
            }
        }
    }

    private static void CheckNode<TIdx>(TIdx idx, IndexContractOracle o, int[] w, IndexLocus locus, int[] alphabet, string context)
        where TIdx : ISuffixTreeNavigation, ISuffixArrayAccess
    {
        string at = $"{context} node {P(w)}";
        Eq(locus, o.LocusOf(w), () => at + " locus");
        IsTrue(idx.IsExplicitNode(locus), () => at + " IsExplicitNode");
        Eq(locus.Length, (long)w.Length, () => at + " Length");
        Eq(locus.Aux, (long)w.Length, () => at + " Aux of a node is its depth");
        CheckCountAndLocate(idx, o, w, locus, at);
        CheckParent(idx, o, w, locus, at);
        CheckSuffixLink(idx, o, w, locus, at);
        CheckEdgeSymbols(idx, o, w, locus, at);
        CheckChildren(idx, o, w, locus, alphabet, at);
        CheckExtensions(idx, o, w, locus, alphabet, at);
    }

    private static void CheckEdge<TIdx>(TIdx idx, IndexContractOracle o, int[] w, IndexLocus node, int[] childString, IndexLocus childLocus, int[] alphabet, string context, CheckOptions options, SplitMix64 rng)
        where TIdx : ISuffixTreeNavigation, ISuffixArrayAccess
    {
        int[] edge = childString[w.Length..];
        int interior = edge[^1] == -1 ? edge.Length - 1 : edge.Length; // symbols of the edge that a pattern can contain
        IndexLocus previous = node;
        var points = new HashSet<int>();
        if (interior <= options.EdgePointCap) for (int k = 1; k <= interior; k++) points.Add(k);
        else
        {
            foreach (int k in new[] { 1, 2, interior / 2, interior - 1, interior }) points.Add(k);
            for (int i = 0; i < 3; i++) points.Add(1 + rng.NextInt(interior));
        }
        int previousK = 0;
        foreach (int k in points.OrderBy(x => x))
        {
            int[] prefix = [.. w, .. edge[..k]];
            string at = $"{context} edge {P(w)} -> {P(childString)} at {k}";
            string pattern = string.Concat(edge[..k].Select(s => (char)s));
            bool ok = idx.TryExtendRight(node, pattern.AsSpan(), out IndexLocus point);
            IsTrue(ok, () => at + " TryExtendRight(pattern)");
            IndexLocus expected = o.LocusOf(prefix);
            Eq(point, expected, () => at + " point");
            if (k == previousK + 1)
            {
                bool stepped = idx.TryExtendRight(previous, edge[k - 1], out IndexLocus viaSymbol);
                IsTrue(stepped, () => at + " TryExtendRight(symbol)");
                Eq(viaSymbol, point, () => at + " symbol step equals pattern jump");
            }
            if (k == edge.Length)
            {
                Eq(point, childLocus, () => at + " the full edge reaches the child");
            }
            else
            {
                IsTrue(!idx.IsExplicitNode(point), () => at + " mid-edge point is no node");
                Eq(point.Length, (long)prefix.Length, () => at + " Length");
                IsTrue(point.Aux > point.Length, () => at + " Aux below the point");
                CheckCountAndLocate(idx, o, prefix, point, at);
                IsTrue(idx.TryParent(point, out IndexLocus parent), () => at + " TryParent");
                Eq(parent, node, () => at + " the parent of an edge point is the top of the edge");
                CheckExtensions(idx, o, prefix, point, alphabet, at);
            }
            previous = point;
            previousK = k;
        }
    }

    private static void CheckCountAndLocate<TIdx>(TIdx idx, IndexContractOracle o, int[] w, IndexLocus locus, string at)
        where TIdx : ISuffixTreeNavigation, ISuffixArrayAccess
    {
        List<int> expected = o.Positions(w).OrderBy(p => o.Isa[p]).ToList(); // suffix-array order, derived from the oracle's ISA
        Eq(idx.Count(locus), (long)expected.Count, () => at + " Count");
        foreach (int page in new[] { 1, 3, 1000 })
        {
            var got = new List<long>();
            var buffer = new long[page];
            long skip = 0;
            while (true)
            {
                int written = idx.Locate(locus, skip, buffer);
                if (written == 0) break;
                got.AddRange(buffer.Take(written));
                skip += written;
            }
            IsTrue(got.SequenceEqual(expected.Select(x => (long)x)), () => $"{at} Locate pages of {page}: got [{string.Join(",", got)}] expected [{string.Join(",", expected)}]");
        }
        if (expected.Count > 2)
        {
            var buffer = new long[2];
            Eq(idx.Locate(locus, 1, buffer), 2, () => at + " Locate(skip 1, 2)");
            IsTrue(buffer.SequenceEqual(expected.Skip(1).Take(2).Select(x => (long)x)), () => at + " Locate skip content");
        }
    }

    private static void CheckParent<TIdx>(TIdx idx, IndexContractOracle o, int[] w, IndexLocus locus, string at)
        where TIdx : ISuffixTreeNavigation, ISuffixArrayAccess
    {
        int[]? expected = o.Parent(w);
        bool got = idx.TryParent(locus, out IndexLocus parent);
        Eq(got, expected is not null, () => at + " TryParent");
        Eq(parent, expected is null ? default : o.LocusOf(expected), () => at + " parent locus");
    }

    private static void CheckSuffixLink<TIdx>(TIdx idx, IndexContractOracle o, int[] w, IndexLocus locus, string at)
        where TIdx : ISuffixTreeNavigation, ISuffixArrayAccess
    {
        int[]? expected = IndexContractOracle.SuffixLink(w);
        bool isInternal = expected is not null && o.Children(w).Count > 0;
        bool got = idx.TrySuffixLink(locus, out IndexLocus link);
        Eq(got, isInternal, () => at + " TrySuffixLink");
        if (isInternal)
        {
            IsTrue(o.IsExplicit(expected!), () => at + " theory: the suffix link target is a node");
            Eq(link, o.LocusOf(expected!), () => at + " suffix link locus");
        }
        else Eq(link, default, () => at + " no link leaves the result default");
    }

    private static void CheckEdgeSymbols<TIdx>(TIdx idx, IndexContractOracle o, int[] w, IndexLocus locus, string at)
        where TIdx : ISuffixTreeNavigation, ISuffixArrayAccess
    {
        int[]? parent = o.Parent(w);
        int parentDepth = parent?.Length ?? 0;
        int edge = w.Length - parentDepth;
        for (int k = 0; k < edge; k++) Eq(idx.EdgeSymbol(locus, k), w[parentDepth + k], () => $"{at} EdgeSymbol({k})");
        if (w.Length == 0 || w[^1] == -1) // one occurrence only (or the empty suffix): the symbols behind the edge are determined as well
            for (int k = edge; k < edge + 3; k++) Eq(idx.EdgeSymbol(locus, k), -1, () => $"{at} EdgeSymbol({k}) behind the terminator");
        Eq(idx.EdgeSymbol(locus, o.N), -1, () => at + " EdgeSymbol(n)");
        Eq(idx.EdgeSymbol(locus, 1L << 40), -1, () => at + " EdgeSymbol(2^40)");
        Eq(idx.EdgeSymbol(locus, long.MaxValue), -1, () => at + " EdgeSymbol(long.MaxValue)");
    }

    private static void CheckChildren<TIdx>(TIdx idx, IndexContractOracle o, int[] w, IndexLocus locus, int[] alphabet, string at)
        where TIdx : ISuffixTreeNavigation, ISuffixArrayAccess
    {
        var kids = o.Children(w);
        int k = kids.Count;
        var exact = new IndexLocus[k + 2];
        Array.Fill(exact, Sentinel);
        Eq(idx.GetChildren(locus, exact), k, () => at + " GetChildren count");
        for (int i = 0; i < k; i++) Eq(exact[i], o.LocusOf(kids[i].Child), () => $"{at} GetChildren[{i}]");
        Eq(exact[k], Sentinel, () => at + " GetChildren wrote behind the last child");
        Eq(exact[k + 1], Sentinel, () => at + " GetChildren wrote far behind the last child");

        var small = new IndexLocus[Math.Max(0, k - 1)];
        Array.Fill(small, Sentinel);
        Eq(idx.GetChildren(locus, small), k, () => at + " GetChildren with a too small span returns the count");
        IsTrue(small.All(x => x == Sentinel), () => at + " GetChildren with a too small span writes nothing");
        Eq(idx.GetChildren(locus, Span<IndexLocus>.Empty), k, () => at + " GetChildren with an empty span");
        var justRight = new IndexLocus[k];
        Eq(idx.GetChildren(locus, justRight), k, () => at + " GetChildren with an exact span");
        for (int i = 0; i < k; i++) Eq(justRight[i], o.LocusOf(kids[i].Child), () => $"{at} GetChildren exact[{i}]");

        var present = kids.Select(x => x.Key).ToHashSet();
        foreach (int key in alphabet.Concat([-1, -2, 65536, int.MaxValue, int.MinValue]))
        {
            if (present.Contains(key)) continue;
            bool found = idx.TryChild(locus, key, out IndexLocus child);
            IsTrue(!found, () => $"{at} TryChild({key}) must not find a child");
            Eq(child, default, () => $"{at} TryChild({key}) leaves the result default");
        }
    }

    private static void CheckExtensions<TIdx>(TIdx idx, IndexContractOracle o, int[] w, IndexLocus locus, int[] alphabet, string at)
        where TIdx : ISuffixTreeNavigation, ISuffixArrayAccess
    {
        foreach (int c in alphabet)
        {
            int[]? right = o.ExtendRight(w, c);
            bool gotRight = idx.TryExtendRight(locus, c, out IndexLocus r);
            Eq(gotRight, right is not null, () => $"{at} TryExtendRight({c})");
            Eq(r, right is null ? default : o.LocusOf(right), () => $"{at} TryExtendRight({c}) locus");
            bool gotPattern = idx.TryExtendRight(locus, new[] { (char)c }, out IndexLocus rp);
            Eq(gotPattern, gotRight, () => $"{at} TryExtendRight(pattern {c})");
            Eq(rp, r, () => $"{at} TryExtendRight(pattern {c}) locus");

            int[]? left = o.ExtendLeft(c, w);
            bool gotLeft = idx.TryExtendLeft(locus, c, out IndexLocus l);
            Eq(gotLeft, left is not null, () => $"{at} TryExtendLeft({c})");
            Eq(l, left is null ? default : o.LocusOf(left), () => $"{at} TryExtendLeft({c}) locus");
        }
        foreach (int bad in new[] { -1, -2, 65536, int.MaxValue, int.MinValue })
        {
            IsTrue(!idx.TryExtendRight(locus, bad, out IndexLocus r), () => $"{at} TryExtendRight({bad})");
            Eq(r, default, () => $"{at} TryExtendRight({bad}) leaves default");
            IsTrue(!idx.TryExtendLeft(locus, bad, out IndexLocus l), () => $"{at} TryExtendLeft({bad})");
            Eq(l, default, () => $"{at} TryExtendLeft({bad}) leaves default");
        }
    }
}
