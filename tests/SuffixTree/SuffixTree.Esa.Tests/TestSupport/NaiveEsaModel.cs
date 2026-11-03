namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>Records the <see cref="ISuffixTreeVisitor"/> events of a traversal as comparable strings.</summary>
public sealed class TraverseRecorder : ISuffixTreeVisitor
{
    /// <summary>The events in call order: <c>N start end leaves children depth</c>, <c>E key</c>, <c>X</c>.</summary>
    public List<string> Events { get; } = new();

    /// <inheritdoc />
    public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth) =>
        Events.Add($"N {startIndex} {endIndex} {leafCount} {childCount} {depth}");

    /// <inheritdoc />
    public void EnterBranch(int key) => Events.Add($"E {key}");

    /// <inheritdoc />
    public void ExitBranch() => Events.Add("X");
}

/// <summary>
/// Test-only naive enhanced-suffix-array tree model for hypothesis H-OCC: the lcp-interval tree of the naive SA/LCP, emitted in
/// the <c>Traverse</c> convention of the Ukkonen tree (children by ascending key with the terminator, key -1, first; a leaf has
/// <c>End = -1</c>; <c>depth</c> is the string depth of the parent; the root is <c>(0, 0, leaves, children, 0)</c>). The edge of
/// an internal node is labelled by the LEFTMOST occurrence of its string (<c>min SA[lb..rb]</c>); the edge of a leaf by its own suffix.
/// The empty text has only the root, as in the Ukkonen tree.
/// </summary>
public static class NaiveEsaModel
{
    /// <summary>Emits the node/branch events of the terminated suffix tree of <paramref name="text"/>.</summary>
    public static void EmitTraverse(string text, ISuffixTreeVisitor visitor)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(visitor);
        int n = text.Length;
        if (n == 0)
        {
            visitor.VisitNode(0, 0, 0, 0, 0);
            return;
        }
        int[] sa = NaiveSuffixArray.Build(text);
        int[] lcp = NaiveLcp.Build(text, sa);
        var groups = Partition(text, sa, 0, n, 0);
        visitor.VisitNode(0, 0, n + 1, groups.Count, 0);
        foreach (var (key, lb, rb) in groups)
            EmitChild(text, sa, lcp, visitor, key, lb, rb, 0);
    }

    private static void EmitChild(string text, int[] sa, int[] lcp, ISuffixTreeVisitor v, int key, int lb, int rb, int parentDepth)
    {
        v.EnterBranch(key);
        if (lb == rb)
        {
            v.VisitNode(sa[lb] + parentDepth, -1, 1, 0, parentDepth);
            v.ExitBranch();
            return;
        }
        int depth = int.MaxValue;
        int occurrence = int.MaxValue;
        for (int r = lb; r <= rb; r++)
        {
            occurrence = Math.Min(occurrence, sa[r]);
            if (r > lb) depth = Math.Min(depth, lcp[r]);
        }
        var kids = Partition(text, sa, lb, rb, depth);
        v.VisitNode(occurrence + parentDepth, occurrence + depth, rb - lb + 1, kids.Count, parentDepth);
        foreach (var (k, l, r) in kids)
            EmitChild(text, sa, lcp, v, k, l, r, depth);
        v.ExitBranch();
    }

    private static List<(int Key, int Lb, int Rb)> Partition(string text, int[] sa, int lb, int rb, int depth)
    {
        var groups = new List<(int, int, int)>();
        int n = text.Length;
        int start = lb;
        int KeyOf(int r) => sa[r] + depth < n ? text[sa[r] + depth] : -1;
        for (int r = lb + 1; r <= rb + 1; r++)
        {
            if (r > rb || KeyOf(r) != KeyOf(start))
            {
                groups.Add((KeyOf(start), start, r - 1));
                start = r;
            }
        }
        return groups;
    }
}
