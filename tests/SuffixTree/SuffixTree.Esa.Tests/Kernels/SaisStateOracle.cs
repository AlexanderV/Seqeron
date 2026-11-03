namespace SuffixTree.Esa.Tests.Kernels;

/// <summary>
/// The intermediate state of one SA-IS recursion level, in the form that every kernel reports and the oracle derives from the definitions:
/// symbol counts and the two bucket arrays of the n + 1 layout (<c>Heads[c]</c> = first index of bucket c, starting at 1; <c>Tails[c]</c> = exclusive end),
/// the LMS positions in text order, the LMS positions in the order of their LMS substrings (after induction round 1), the number of names and the reduced string.
/// </summary>
internal sealed record SaisLevelState(
    int Length,
    int Sigma,
    int[] Counts,
    int[] Heads,
    int[] Tails,
    int[] LmsPositions,
    int[] LmsSorted,
    int NameCount,
    int[] Reduced);

/// <summary>
/// Test oracle for the intermediate states of SA-IS (review MAJ-1 / section 6 T1): everything is computed from the definitions (types, LMS
/// positions, LMS substrings compared as sequences of (symbol, type), dense ranks as names) and shares no code with <see cref="SaisReference"/>.
/// The recursion follows the algorithm: a level is reduced again iff it has fewer names than LMS positions.
/// </summary>
internal static class SaisStateOracle
{
    /// <summary>The states of all levels in recursion order for the symbol sequence <paramref name="text"/> over <c>[0, sigma)</c>.</summary>
    internal static List<SaisLevelState> Compute(int[] text, int sigma)
    {
        var levels = new List<SaisLevelState>();
        int[] current = text;
        int currentSigma = sigma;
        while (true)
        {
            SaisLevelState state = Level(current, currentSigma);
            levels.Add(state);
            if (state.NameCount >= state.LmsPositions.Length) return levels; // all names distinct (or no LMS position): direct placement, no recursion
            current = state.Reduced;
            currentSigma = state.NameCount;
        }
    }

    /// <summary>The name (rank among the distinct LMS substrings) of every LMS position of a level, indexed by text position (-1 for non-LMS positions).</summary>
    internal static int[] NameOfPosition(SaisLevelState level)
    {
        var names = new int[level.Length + 1];
        Array.Fill(names, -1);
        for (int i = 0; i < level.LmsPositions.Length; i++) names[level.LmsPositions[i]] = level.Reduced[i];
        return names;
    }

    private static SaisLevelState Level(int[] t, int sigma)
    {
        int n = t.Length;
        var counts = new int[sigma];
        foreach (int c in t) counts[c]++;
        var heads = new int[sigma];
        var tails = new int[sigma];
        int sum = 1;
        for (int c = 0; c < sigma; c++)
        {
            heads[c] = sum;
            sum += counts[c];
            tails[c] = sum;
        }

        var isS = new bool[n + 1];
        isS[n] = true;
        for (int i = n - 2; i >= 0; i--) isS[i] = t[i] < t[i + 1] || (t[i] == t[i + 1] && isS[i + 1]);
        var lms = new List<int>();
        for (int i = 1; i < n; i++)
            if (isS[i] && !isS[i - 1]) lms.Add(i);

        // LMS substring of p: positions p..q inclusive where q is the next LMS position or the sentinel n; elements are (symbol, type), sentinel = -1
        var keys = new Dictionary<int, int[]>();
        for (int k = 0; k < lms.Count; k++)
        {
            int p = lms[k];
            int q = k + 1 < lms.Count ? lms[k + 1] : n;
            var key = new List<int>();
            for (int x = p; x <= q; x++) key.Add(x == n ? -1 : (t[x] * 2) + (isS[x] ? 1 : 0));
            keys[p] = [.. key];
        }

        int[] sorted = [.. lms.OrderBy(p => keys[p], Comparer<int[]>.Create(CompareKeys)).ThenBy(p => p)];
        var nameOf = new Dictionary<int, int>();
        int name = -1;
        for (int i = 0; i < sorted.Length; i++)
        {
            if (i == 0 || CompareKeys(keys[sorted[i - 1]], keys[sorted[i]]) != 0) name++;
            nameOf[sorted[i]] = name;
        }
        int[] reduced = [.. lms.Select(p => nameOf[p])];
        return new SaisLevelState(n, sigma, counts, heads, tails, [.. lms], sorted, name + 1, reduced);
    }

    private static int CompareKeys(int[] a, int[] b)
    {
        int k = Math.Min(a.Length, b.Length);
        for (int i = 0; i < k; i++)
            if (a[i] != b[i]) return a[i].CompareTo(b[i]);
        return a.Length.CompareTo(b.Length);
    }
}
