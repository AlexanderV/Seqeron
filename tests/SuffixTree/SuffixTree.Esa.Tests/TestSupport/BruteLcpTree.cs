namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>One lcp-interval of the definitional tree (test oracle L0): <c>[Lb..Rb]</c> with value <c>Ell</c>, its l-indices in ascending order, and whether it is the last child of its parent.</summary>
public sealed record BruteInterval(int Ell, int Lb, int Rb, IReadOnlyList<int> LIndices, bool LastChild);

/// <summary>
/// The lcp-interval tree built straight from the definition by recursion on the LCP array, without a stack scan and without a child table:
/// the interval <c>[lb..rb]</c> has the value <c>ell = min LCP[lb+1..rb]</c>, its l-indices are the positions of <c>(lb, rb]</c> where the LCP equals
/// <c>ell</c>, and its children start at <c>lb</c> and at every l-index. An explicit work list replaces the recursion (a^n has depth n).
/// </summary>
public static class BruteLcpTree
{
    /// <summary>Every internal interval (the root included) of the terminated tree for the n + 1 entry LCP array <paramref name="lcp"/>; empty for n = 0.</summary>
    public static List<BruteInterval> Intervals(int[] lcp)
    {
        ArgumentNullException.ThrowIfNull(lcp);
        int n = lcp.Length - 1;
        var result = new List<BruteInterval>();
        if (n < 1) return result;
        var work = new Stack<(int Lb, int Rb, bool Last)>();
        work.Push((0, n, true));
        while (work.Count > 0)
        {
            var (lb, rb, last) = work.Pop();
            int ell = int.MaxValue;
            for (int k = lb + 1; k <= rb; k++) ell = Math.Min(ell, lcp[k]);
            var lIndices = new List<int>();
            for (int k = lb + 1; k <= rb; k++)
                if (lcp[k] == ell) lIndices.Add(k);
            result.Add(new BruteInterval(ell, lb, rb, lIndices, last));
            int start = lb;
            for (int t = 0; t <= lIndices.Count; t++)
            {
                int end = t < lIndices.Count ? lIndices[t] - 1 : rb;
                if (end > start) work.Push((start, end, t == lIndices.Count));
                start = end + 1;
            }
        }
        return result;
    }
}
