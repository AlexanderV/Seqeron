using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// The checks of the long navigation engine over a <see cref="NavRig"/> (work package S2.5), each against the closed form and nothing of the engine: the suffix link of the internal node of a factor is the node of the factor without its
/// first symbol, the parent of a node is the node of the prefix of its parent's depth, the Lcp of two ranks is the range minimum of the closed form, the Lca is the node of the common prefix, and the children of a node partition its
/// interval in the order of their first symbols. The default set runs them at n of a few thousand over every family; the <c>Large</c> category at 2^32 + 10 and 2^33.
/// </summary>
internal static class NavChecks
{
    /// <summary>The suffix link of every internal node among the factors of the rig; returns how many nodes were checked.</summary>
    internal static int SuffixLinks(NavRig rig, string label)
    {
        int checkedNodes = 0;
        foreach ((long start, long length) in rig.Factors())
        {
            EsaNode64? node = rig.NodeOfFactor(start, length);
            if (node is not { } v) continue;
            EsaNode64 link = rig.Nav.SuffixLink(v);
            if (length == 1)
            {
                Assert.That(link, Is.EqualTo(rig.Nav.Root), $"{label} the link of a node of depth 1 is the root ({start}, {length})");
            }
            else
            {
                RankRange expected = rig.Form.RangeOfFactor(start + 1, length - 1);
                Assert.That((link.Lb, link.Rb, link.Depth), Is.EqualTo((expected.First, expected.Last, length - 1)), $"{label} link of the factor ({start}, {length})");
                EsaNode64 full = rig.NodeOf(expected.First, expected.Last, length - 1);
                Assert.That(link, Is.EqualTo(full), $"{label} the link carries its parent depth, last-child flag and position ({start}, {length})");
            }
            checkedNodes++;
        }
        return checkedNodes;
    }

    /// <summary>The parent of every internal node among the factors.</summary>
    internal static int Parents(NavRig rig, string label)
    {
        int checkedNodes = 0;
        foreach ((long start, long length) in rig.Factors())
        {
            EsaNode64? node = rig.NodeOfFactor(start, length);
            if (node is not { } v) continue;
            EsaNode64 parent = rig.Nav.Parent(v);
            if (v.ParentDepth == 0)
            {
                Assert.That(parent, Is.EqualTo(rig.Nav.Root), $"{label} the parent of ({start}, {length}) is the root");
            }
            else
            {
                RankRange expected = rig.Form.RangeOfFactor(v.Pos, v.ParentDepth);
                Assert.That((parent.Lb, parent.Rb, parent.Depth), Is.EqualTo((expected.First, expected.Last, v.ParentDepth)), $"{label} parent of the factor ({start}, {length})");
                Assert.That(parent, Is.EqualTo(rig.NodeOf(expected.First, expected.Last, v.ParentDepth)), $"{label} the parent carries its fields ({start}, {length})");
            }
            checkedNodes++;
        }
        return checkedNodes;
    }

    /// <summary>The Lcp and the Lca of pairs of ranks: the ends, the middle, the boundary 2^31 and seeded random ranks.</summary>
    internal static void LcpAndLca(NavRig rig, string label, int randomPairs)
    {
        long n = rig.Form.Length;
        var ranks = new SortedSet<long> { 0, 1, 2, 3, n - 3, n - 2, n - 1, n, n / 2, n / 3, (1L << 31) - 1, 1L << 31, (1L << 31) + 1 };
        var rng = new SplitMix64(0x5EED25 + (ulong)n);
        for (int i = 0; i < randomPairs; i++) ranks.Add((long)(rng.Next() % (ulong)(n + 1)));
        long[] list = [.. ranks.Where(r => r >= 0 && r <= n)];
        foreach (long a in list)
        {
            foreach (long b in list.Where(x => x >= a).Take(24).Concat([a]).Distinct())
            {
                long lcp = a == b ? n - rig.Form.SuffixAt(a) : rig.Form.LcpMin(Math.Min(a, b) + 1, Math.Max(a, b));
                Assert.That(rig.Nav.Lcp(a, b), Is.EqualTo(lcp), $"{label} Lcp({a}, {b})");
                Assert.That(rig.Nav.Lcp(b, a), Is.EqualTo(lcp), $"{label} Lcp({b}, {a})");
                (long lb, long rb, long depth) = rig.Nav.Lca(a, b);
                if (n == 0) continue;
                if (a == b)
                {
                    Assert.That((lb, rb, depth), Is.EqualTo((a, a, n - rig.Form.SuffixAt(a) + 1)), $"{label} Lca({a}, {a}) is the leaf");
                }
                else if (lcp == 0)
                {
                    Assert.That((lb, rb, depth), Is.EqualTo((0L, n, 0L)), $"{label} Lca({a}, {b}) is the root");
                }
                else
                {
                    RankRange expected = rig.Form.RangeOfFactor(rig.Form.SuffixAt(a), lcp);
                    Assert.That((lb, rb, depth), Is.EqualTo((expected.First, expected.Last, lcp)), $"{label} Lca({a}, {b})");
                }
            }
        }
    }

    /// <summary>The children of every internal node among the factors partition its interval, in the order of their first symbols, each an internal node or a leaf of the right depth.</summary>
    internal static int Children(NavRig rig, string label)
    {
        int checkedNodes = 0;
        long n = rig.Form.Length;
        foreach ((long start, long length) in rig.Factors())
        {
            EsaNode64? node = rig.NodeOfFactor(start, length);
            if (node is not { } v) continue;
            var list = new List<EsaNode64>();
            var sink = new ListSink(list);
            int count = rig.Nav.Children(v, ref sink);
            Assert.That(count, Is.EqualTo(list.Count));
            Assert.That(count, Is.GreaterThanOrEqualTo(2), $"{label} an internal node has two children or more ({start}, {length})");
            long next = v.Lb;
            int previousKey = int.MinValue;
            foreach (EsaNode64 child in list)
            {
                Assert.That(child.Lb, Is.EqualTo(next), $"{label} the children follow each other ({start}, {length})");
                Assert.That(child.ParentDepth, Is.EqualTo(v.Depth));
                if (child.Lb == child.Rb) Assert.That(child.Depth, Is.EqualTo(n - child.Pos + 1), $"{label} a leaf has the depth of its suffix and the terminator");
                else Assert.That(rig.IsInternalNode(child.Lb, child.Rb, child.Depth), Is.True, $"{label} an interval child is a node ({start}, {length})");
                Assert.That(child, Is.EqualTo(rig.NodeOf(child.Lb, child.Rb, child.Depth)), $"{label} the child carries its fields");
                long at = child.Pos + v.Depth;
                int key = at < n ? rig.Form.CharAt(at) : -1;
                Assert.That(key, Is.GreaterThan(previousKey), $"{label} the keys ascend ({start}, {length})");
                Assert.That(rig.Nav.TryChild(v, key, out EsaNode64 byKey), Is.True);
                Assert.That(byKey, Is.EqualTo(child), $"{label} TryChild by the key of the child");
                previousKey = key;
                next = child.Rb + 1;
            }
            Assert.That(next, Is.EqualTo(v.Rb + 1), $"{label} the children end at the end of the node ({start}, {length})");
            Assert.That(rig.Nav.TryChild(v, 'z', out _), Is.False);
            checkedNodes++;
        }
        return checkedNodes;
    }
}
