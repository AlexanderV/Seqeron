using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Heavy-tier properties for the REP-DIRECT-001 enumeration variants added by review batch B04 (audit WP2):
/// reverse-complement maximal pairs (MUMmer repeat-match / Vmatch -p), maximal k-mismatch repeats
/// (REPuter / Vmatch -h k, both maximality readings) and supermaximal repeats (Gusfield §7.12.1 / Vmatch -supermax),
/// each equal to an independent brute force of the published definition.
/// Test Unit: REP-DIRECT-001.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Repeats")]
public class RepDirectVariantsProperties
{
    private static Gen<string> Dna(string alphabets, int minLen, int maxLen) =>
        from alpha in Gen.Elements(alphabets.Split('|'))
        from n in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements(alpha.ToCharArray()).ArrayOf(n)
        select new string(chars);

    private static bool Acgt(char c) => "ACGTacgt".Contains(c);

    private static bool Pairs(char a, char b) =>
        Acgt(a) && char.ToUpperInvariant(SequenceExtensions.GetComplementBase(a)) == char.ToUpperInvariant(b);

    private static bool Same(char a, char b) => Acgt(a) && char.ToUpperInvariant(a) == char.ToUpperInvariant(b);

    internal static List<(int, int, int)> BruteReverseComplement(string s, int minL, int maxL, int minSp)
    {
        var o = new List<(int, int, int)>();
        int n = s.Length;
        for (int i = 0; i < n; i++)
            for (int k = i; k < n; k++)
                for (int len = 1; i + len <= n && k + len <= n; len++)
                {
                    bool ok = true;
                    for (int t = 0; t < len && ok; t++) ok = Pairs(s[i + t], s[k + len - 1 - t]);
                    if (!ok) continue;
                    bool outward = i > 0 && k + len < n && Pairs(s[i - 1], s[k + len]);
                    bool inward = i + len < n && k > 0 && Pairs(s[i + len], s[k - 1]);
                    if (!outward && !inward && len >= minL && len <= maxL && (long)k - i - len >= minSp)
                        o.Add((i, k, len));
                }
        return o.OrderBy(t => t.Item1).ThenBy(t => t.Item2).ThenBy(t => t.Item3).ToList();
    }

    internal static List<(int, int, int, int)> BruteMismatch(string s, int minL, int k, bool literal)
    {
        var all = new List<(int, int, int, int)>();
        int n = s.Length;
        for (int d = 1; d < n; d++)
        {
            int m = n - d;
            for (int st = 0; st < m; st++)
            {
                int c = 0;
                for (int e = st; e < m; e++)
                {
                    c += Same(s[e], s[e + d]) ? 0 : 1;
                    if (c > k) break;
                    bool lext = st > 0 && c + (Same(s[st - 1], s[st - 1 + d]) ? 0 : 1) <= k;
                    bool rext = e + 1 < m && c + (Same(s[e + 1], s[e + 1 + d]) ? 0 : 1) <= k;
                    if (!lext && !rext && e - st + 1 >= minL) all.Add((st, st + d, e - st + 1, c));
                }
            }
        }
        if (literal)
        {
            all = all.Where(a => !all.Any(b => b != a && b.Item1 <= a.Item1 && a.Item1 + a.Item3 <= b.Item1 + b.Item3
                                               && b.Item2 <= a.Item2 && a.Item2 + a.Item3 <= b.Item2 + b.Item3)).ToList();
        }
        return all.OrderBy(t => t.Item1).ThenBy(t => t.Item2).ThenBy(t => t.Item3).ToList();
    }

    internal static List<(int, string)> BruteSupermaximal(string s, int minL)
    {
        string u = s.ToUpperInvariant();
        int n = u.Length;
        var reps = new HashSet<string>();
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
            {
                if (i > 0 && Same(u[i - 1], u[j - 1])) continue;
                int len = 0;
                while (j + len < n && Same(u[i + len], u[j + len])) len++;
                if (len > 0) reps.Add(u.Substring(i, len));
            }
        var result = new List<(int, string)>();
        foreach (var r in reps.Where(r => r.Length >= minL && !reps.Any(o => o != r && o.Contains(r, StringComparison.Ordinal))))
        {
            var pos = Enumerable.Range(0, n - r.Length + 1).Where(p => string.CompareOrdinal(u, p, r, 0, r.Length) == 0);
            result.Add((r.Length, string.Join(",", pos)));
        }
        return result.OrderBy(t => t.Item2).ThenBy(t => t.Item1).ToList();
    }

    /// <summary>Reverse-complement maximal pairs equal the brute force of the Vmatch/repeat-match definition.</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property ReverseComplement_EqualsBruteForce()
    {
        var gen = from s in Dna("ACGT|AT|GC|ACGTN|acgT|ACGU", 0, 40)
                  from minL in Gen.Choose(2, 6)
                  from maxExtra in Gen.Elements(0, 3, int.MaxValue - 10)
                  from minSp in Gen.Elements(int.MinValue, -4, 0, 1, 3)
                  select (s, minL, maxL: maxExtra > 1000 ? int.MaxValue : minL + maxExtra, minSp);
        return Prop.ForAll(gen.ToArbitrary(), c =>
        {
            var expected = BruteReverseComplement(c.s, c.minL, c.maxL, c.minSp);
            var actual = RepeatFinder.FindReverseComplementRepeats(c.s, c.minL, c.maxL, c.minSp)
                .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();
            return actual.SequenceEqual(expected).Label($"{c}");
        });
    }

    /// <summary>Maximal k-mismatch repeats (per-diagonal and literal-containment readings) equal the brute force.</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property KMismatch_EqualsBruteForce()
    {
        var gen = from s in Dna("ACGT|AC|AT|ACGTN|acgt", 0, 40)
                  from k in Gen.Choose(0, 3)
                  from extra in Gen.Choose(1, 8)
                  from literal in Gen.Elements(true, false)
                  select (s, k, minL: Math.Max(2, k + extra), literal);
        return Prop.ForAll(gen.ToArbitrary(), c =>
        {
            var expected = BruteMismatch(c.s, c.minL, c.k, c.literal);
            var actual = RepeatFinder.FindApproximateDirectRepeats(c.s, c.minL, c.k, int.MaxValue, int.MinValue, c.literal)
                .Select(r => (r.FirstPosition, r.SecondPosition, r.Length, r.Mismatches)).ToList();
            return actual.SequenceEqual(expected).Label($"{c}");
        });
    }

    /// <summary>Supermaximal repeats equal the brute force "maximal repeats contained in no other maximal repeat".</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property Supermaximal_EqualsBruteForce()
    {
        var gen = from s in Dna("ACGT|AC|AT|ACGTN|acgt|NNA", 0, 50)
                  from minL in Gen.Choose(1, 6)
                  select (s, minL);
        return Prop.ForAll(gen.ToArbitrary(), c =>
        {
            var expected = BruteSupermaximal(c.s, c.minL);
            var actual = RepeatFinder.FindSupermaximalRepeats(c.s, c.minL)
                .Select(r => (r.Length, string.Join(",", r.Positions))).OrderBy(t => t.Item2).ThenBy(t => t.Item1).ToList();
            return actual.SequenceEqual(expected).Label($"{c}");
        });
    }
}
