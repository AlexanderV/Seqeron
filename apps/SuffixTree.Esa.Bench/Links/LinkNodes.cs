namespace SuffixTree.Esa.Bench.Links;

/// <summary>
/// All internal nodes (lcp-intervals of value &gt;= 1; the root is not a link source) of an index, with string depth and node depth (number of edges from the root).
/// Built by one stack scan over the LCP array (an independent re-implementation of <see cref="LcpIntervalScan"/>, checked against the definitional interval tree).
/// </summary>
internal sealed class LinkNodeSet
{
    private LinkNodeSet(int[] lb, int[] rb, int[] depth, int[] nodeDepth)
    {
        Lb = lb;
        Rb = rb;
        Depth = depth;
        NodeDepth = nodeDepth;
    }

    /// <summary>Number of internal nodes.</summary>
    internal int Count => Lb.Length;

    /// <summary>Left boundaries.</summary>
    internal int[] Lb { get; }

    /// <summary>Right boundaries.</summary>
    internal int[] Rb { get; }

    /// <summary>String depths (the lcp-interval values).</summary>
    internal int[] Depth { get; }

    /// <summary>Node depths: 1 for a child of the root.</summary>
    internal int[] NodeDepth { get; }

    /// <summary>Enumerates the internal nodes by one stack scan over the LCP array.</summary>
    internal static LinkNodeSet Enumerate(CompactLcp lcp)
    {
        ArgumentNullException.ThrowIfNull(lcp);
        int n = lcp.Length - 1;
        int capacity = Math.Max(16, Math.Min(n + 1, 1 << 16));
        int[] eLb = new int[capacity], eRb = new int[capacity], eEll = new int[capacity], eParentPush = new int[capacity];
        int[] emissionOfPush = new int[capacity];
        int count = 0, pushes = 0;

        int stackCapacity = 64;
        int[] stEll = new int[stackCapacity], stLb = new int[stackCapacity], stPid = new int[stackCapacity];
        stEll[0] = 0;
        stLb[0] = 0;
        stPid[0] = pushes++;
        int sp = 1;
        var cursor = lcp.CreateCursor();

        for (int i = 1; i <= n + 1; i++)
        {
            int cur = i <= n ? cursor.Get(i) : -1;
            int left = i - 1;
            while (sp > 0 && stEll[sp - 1] > cur)
            {
                sp--;
                if (count == eLb.Length)
                {
                    int grown = eLb.Length * 2;
                    Array.Resize(ref eLb, grown);
                    Array.Resize(ref eRb, grown);
                    Array.Resize(ref eEll, grown);
                    Array.Resize(ref eParentPush, grown);
                }
                if (stPid[sp] >= emissionOfPush.Length) Array.Resize(ref emissionOfPush, Math.Max(emissionOfPush.Length * 2, stPid[sp] + 1));
                eLb[count] = stLb[sp];
                eRb[count] = i - 1;
                eEll[count] = stEll[sp];
                emissionOfPush[stPid[sp]] = count;
                // the parent is the interval below on the stack when its value is not below cur (it is popped later or stays), else the interval pushed for cur right after this loop
                int parentPush = -1;
                if (sp > 0) parentPush = stEll[sp - 1] >= cur ? stPid[sp - 1] : pushes;
                eParentPush[count] = parentPush;
                count++;
                left = stLb[sp];
            }
            if (i > n) break;

            if (stEll[sp - 1] < cur)
            {
                if (sp == stackCapacity)
                {
                    stackCapacity *= 2;
                    Array.Resize(ref stEll, stackCapacity);
                    Array.Resize(ref stLb, stackCapacity);
                    Array.Resize(ref stPid, stackCapacity);
                }
                stEll[sp] = cur;
                stLb[sp] = left;
                stPid[sp] = pushes++;
                sp++;
            }
        }

        // post-order: a parent is emitted after its children, so node depths are final when read in reverse emission order
        var nodeDepthAll = new int[count];
        for (int e = count - 1; e >= 0; e--)
            nodeDepthAll[e] = eParentPush[e] < 0 ? 0 : nodeDepthAll[emissionOfPush[eParentPush[e]]] + 1;

        int internalCount = count - 1;      // the root is emitted last
        if (internalCount < 0) internalCount = 0;
        var lb = new int[internalCount];
        var rb = new int[internalCount];
        var depth = new int[internalCount];
        var nodeDepth = new int[internalCount];
        Array.Copy(eLb, lb, internalCount);
        Array.Copy(eRb, rb, internalCount);
        Array.Copy(eEll, depth, internalCount);
        Array.Copy(nodeDepthAll, nodeDepth, internalCount);
        return new LinkNodeSet(lb, rb, depth, nodeDepth);
    }

    /// <summary>The handle of node <paramref name="i"/>.</summary>
    internal LinkNode Node(int i, int[] sa) => new(Lb[i], Rb[i], Depth[i], sa[Lb[i]]);

    /// <summary>A deterministic random sample of <paramref name="count"/> node indices without repetition (all of them, shuffled, when fewer exist).</summary>
    internal int[] Sample(int count, ulong seed)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(count);
        int take = Math.Min(count, Count);
        var all = new int[Count];
        for (int i = 0; i < all.Length; i++) all[i] = i;
        ulong state = seed;
        for (int i = 0; i < take; i++)
        {
            int j = i + (int)(Next(ref state) % (ulong)(all.Length - i));
            (all[i], all[j]) = (all[j], all[i]);
        }
        return all.AsSpan(0, take).ToArray();
    }

    private static ulong Next(ref ulong state)
    {
        unchecked
        {
            state += 0x9E3779B97F4A7C15UL;
            ulong z = state;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
            return z ^ (z >> 31);
        }
    }

    /// <summary>The indices of the <paramref name="count"/> largest keys, largest first (ties: smaller index first).</summary>
    internal static int[] Largest(ReadOnlySpan<int> keys, int count)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(count);
        count = Math.Min(count, keys.Length);
        if (count == 0) return [];
        // min-heap of the best `count` so far; the top is the worst of them (smallest key, and among equal keys the largest index)
        var heap = new PriorityQueue<int, (int Key, int NegativeIndex)>();
        for (int i = 0; i < keys.Length; i++)
        {
            var priority = (keys[i], -i);
            if (heap.Count < count)
            {
                heap.Enqueue(i, priority);
            }
            else if (heap.TryPeek(out _, out var worst) && priority.CompareTo(worst) > 0)
            {
                heap.Dequeue();
                heap.Enqueue(i, priority);
            }
        }
        var result = new int[heap.Count];
        for (int k = result.Length - 1; k >= 0; k--) result[k] = heap.Dequeue();
        return result;
    }
}
