namespace SuffixTree.Esa;

internal sealed partial class EsaIndex
{
    /// <summary>One open internal node of the iterative walk: the next child starts at <c>A</c>, the one after it at <c>NextBoundary</c> (-1: the next child is the last).</summary>
    private struct TraverseFrame
    {
        internal int A;
        internal int NextBoundary;
        internal int Rb;
        internal int Ell;
    }

    /// <summary>
    /// Walks the lcp-interval tree top-down in suffix-array order and emits the events of the Ukkonen tree
    /// (<c>SuffixTree.Diagnostics.cs</c>, <c>Traverse</c>): <c>VisitNode(start, end, leafCount, childCount, parentDepth)</c> for the root, then for every child in ascending
    /// key order (the terminator, key -1, first) <c>EnterBranch(key)</c>, the child's node and subtree, <c>ExitBranch()</c>. Positions: a leaf has <c>start = SA[i] + parentDepth</c>
    /// and <c>end = -1</c> (boundless); an internal node has <c>start = occ + parentDepth</c> and <c>end = occ + depth</c> with <c>occ = min SA[lb..rb]</c> (hypothesis H-OCC,
    /// <see cref="LeftmostOccurrence"/>); the root is <c>(0, 0, n + 1, children, 0)</c>, and the empty text has the single event <c>(0, 0, 0, 0, 0)</c>. The walk is iterative
    /// (a^n is a chain of n internal nodes) and needs a transient 4 bytes per rank for the occurrences plus one frame per open node.
    /// </summary>
    internal void Traverse(ISuffixTreeVisitor visitor)
    {
        ArgumentNullException.ThrowIfNull(visitor);
        int n = Length;
        if (n == 0)
        {
            visitor.VisitNode(0, 0, 0, 0, 0);
            return;
        }

        SaView sa = Sa;
        EsaText text = TextView;
        bool mapped = IsMapped;
        int[] occurrence = LeftmostOccurrence.Compute(sa, Lcp);
        ulong[]? seen = mapped ? new ulong[(n >> 6) + 1] : null; // one bit per text position: the leaves of a suffix array that is a permutation are distinct (a repeated one would make a visitor that links leaves by position go round in a circle)
        int firstOfRoot = Child.FirstLIndex(0, n, true);
        if (mapped) CheckFirst(firstOfRoot, 0, n);
        visitor.VisitNode(0, 0, n + 1, CountChildren(firstOfRoot, n), 0);

        var frames = new List<TraverseFrame>(64) { new() { A = 0, NextBoundary = firstOfRoot, Rb = n, Ell = 0 } };
        while (frames.Count > 0)
        {
            int top = frames.Count - 1;
            TraverseFrame f = frames[top];
            if (f.A > f.Rb)
            {
                frames.RemoveAt(top);
                if (top > 0) visitor.ExitBranch();
                continue;
            }

            bool lastChild = f.NextBoundary < 0;
            int a = f.A;
            int b = lastChild ? f.Rb : f.NextBoundary - 1;
            int next = lastChild ? -1 : Child.NextLIndex(f.NextBoundary, f.Rb);
            if (mapped && (b < a || b > f.Rb || (next >= 0 && (next <= f.NextBoundary || next > f.Rb)))) throw DamagedChildTable();
            frames[top] = new TraverseFrame { A = b + 1, NextBoundary = next, Rb = f.Rb, Ell = f.Ell };

            int p = sa[a] + f.Ell;
            visitor.EnterBranch((uint)p < (uint)n ? text[p] : -1);
            if (a == b)
            {
                if (seen is not null) MarkLeaf(seen, sa[a], n);
                visitor.VisitNode(p, -1, 1, 0, f.Ell);
                visitor.ExitBranch();
                continue;
            }

            int first = Child.FirstLIndex(a, b, lastChild);
            if (mapped) CheckFirst(first, a, b);
            int ell = Lcp[first];
            if (mapped && (ell <= f.Ell || ell > n)) throw DamagedChildTable();
            int occ = occurrence[first];
            visitor.VisitNode(occ + f.Ell, occ + ell, b - a + 1, CountChildren(first, b), f.Ell);
            frames.Add(new TraverseFrame { A = a, NextBoundary = first, Rb = b, Ell = ell });
        }
    }

    /// <summary>Number of children of an internal interval ending at <paramref name="rb"/> whose first l-index is <paramref name="first"/>: the l-indices plus the first child.</summary>
    private int CountChildren(int first, int rb)
    {
        int count = 1;
        for (int k = first; k >= 0;)
        {
            count++;
            int next = Child.NextLIndex(k, rb);
            if (IsMapped && next >= 0 && (next <= k || next > rb)) throw DamagedChildTable();
            k = next;
        }
        return count;
    }

    /// <summary>The first l-index of an interval that spans <c>[lb..rb]</c> must lie in <c>(lb, rb]</c>.</summary>
    private static void CheckFirst(int first, int lb, int rb)
    {
        if (first <= lb || first > rb) throw DamagedChildTable();
    }

    /// <summary>Records the text position of a leaf; a position outside the text, or one that was reported before, means that the suffix array of the file is not a permutation.</summary>
    private static void MarkLeaf(ulong[] seen, int position, int n)
    {
        if ((uint)position > (uint)n) throw new InvalidDataException("the suffix array holds a position outside the text: the file is damaged");
        ulong bit = 1UL << (position & 63);
        ref ulong word = ref seen[position >> 6];
        if ((word & bit) != 0) throw new InvalidDataException("the suffix array names a text position twice: the file is damaged");
        word |= bit;
    }

    private static InvalidDataException DamagedChildTable() => new("the child table or the LCP array is inconsistent: the file is damaged");
}
