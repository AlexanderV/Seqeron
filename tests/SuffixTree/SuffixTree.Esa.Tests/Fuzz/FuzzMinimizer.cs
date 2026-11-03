namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// Failure minimizer: shrinks a failing text while the failure persists, by deleting chunks of halving size (down to single code units), merging
/// two symbols into one, and finally renaming the symbols to 'a', 'b', ... when that keeps the failure. The result is 1-minimal with respect to
/// the deletion of a single code unit (unless the check budget ran out).
/// </summary>
internal static class FuzzMinimizer
{
    /// <summary>Shrinks <paramref name="text"/>; <paramref name="fails"/> must be deterministic and true for <paramref name="text"/>.</summary>
    internal static string Minimize(string text, Func<string, bool> fails, int maxChecks = 50_000)
    {
        int checks = 0;
        bool Fails(string s)
        {
            checks++;
            return fails(s);
        }

        string cur = text;
        bool changed = true;
        while (changed && checks < maxChecks)
        {
            changed = false;

            for (int chunk = Math.Max(1, cur.Length / 2); chunk >= 1 && checks < maxChecks; chunk /= 2)
            {
                int i = 0;
                while (i < cur.Length && checks < maxChecks)
                {
                    int len = Math.Min(chunk, cur.Length - i);
                    string candidate = cur.Remove(i, len);
                    if (Fails(candidate))
                    {
                        cur = candidate;
                        changed = true;
                    }
                    else
                    {
                        i += len;
                    }
                }
            }

            char[] symbols = cur.Distinct().Order().ToArray();
            for (int x = symbols.Length - 1; x >= 1 && checks < maxChecks; x--)
            {
                for (int y = 0; y < symbols.Length; y++)
                {
                    if (y == x || !cur.Contains(symbols[x])) continue;
                    string candidate = cur.Replace(symbols[x], symbols[y]);
                    if (Fails(candidate))
                    {
                        cur = candidate;
                        changed = true;
                        break;
                    }
                }
            }
        }

        string canonical = Canonical(cur);
        return canonical != cur && Fails(canonical) ? canonical : cur;
    }

    /// <summary>Order-preserving renaming of the distinct symbols to 'a', 'b', ... (unchanged when there are more than 26).</summary>
    private static string Canonical(string text)
    {
        char[] symbols = text.Distinct().Order().ToArray();
        if (symbols.Length > 26) return text;
        var c = new char[text.Length];
        for (int i = 0; i < c.Length; i++) c[i] = (char)('a' + Array.BinarySearch(symbols, text[i]));
        return new string(c);
    }
}
