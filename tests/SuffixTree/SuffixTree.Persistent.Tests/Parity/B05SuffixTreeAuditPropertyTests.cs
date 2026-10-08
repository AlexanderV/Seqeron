namespace SuffixTree.Persistent.Tests.Parity;

/// <summary>
/// Heavy-tier property / metamorphic tests for the review 2026-09 B05 audit additions on SUFFIXTREE-CORE
/// (docs/Validation/review-2026-09/B05.md F26: maximal repeated pairs, all distinct LCS, k-common substrings), run on
/// both the in-memory <see cref="global::SuffixTree.SuffixTree"/> and the persistent tree:
///   • MRP-BRUTE   FindMaximalRepeatedPairs == the Gusfield §7.12 definition by brute force (i &lt; j, equal,
///                 right-maximal, left-maximal; unique symbols match nothing), ascending (first, second).
///   • MRP-MIN     pairs(minLength + d) == pairs(minLength) filtered by Length ≥ minLength + d.
///   • MRP-RELABEL a bijective relabelling of the alphabet leaves the pairs unchanged.
///   • MRP-UNIQUE  a unique-symbol predicate == replacing each such occurrence by a distinct fresh symbol (documented).
///   • LCS-SWAP    the set of distinct LCS strings of (text, other) == that of (other, text); every entry has the
///                 brute-force LCS length and its positions are all occurrences in each string.
///   • KCS-MONO    l(q) is non-increasing in q, l(1) = longest text; permuting the texts changes neither l nor the
///                 per-q substring sets; every reported substring has length l(q) and occurs in ≥ q texts, and no
///                 longer substring does (brute force).
/// Fixed seeds, bounded sizes.
/// </summary>
[TestFixture]
[Category("Parity")]
public class B05SuffixTreeAuditPropertyTests
{
    private static string RandomString(Random rng, int length, string alphabet)
    {
        var c = new char[length];
        for (int i = 0; i < length; i++)
            c[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(c);
    }

    private static List<(int, int, int)> BruteForcePairs(string t, int minLength, Func<char, bool>? unique)
    {
        bool Eq(int x, int y) => t[x] == t[y] && (unique is null || !unique(t[x]));
        var result = new List<(int, int, int)>();
        for (int i = 0; i < t.Length; i++)
            for (int j = i + 1; j < t.Length; j++)
            {
                if (i > 0 && Eq(i - 1, j - 1)) continue; // not left-maximal
                int l = 0;
                while (j + l < t.Length && Eq(i + l, j + l)) l++;
                if (l >= minLength) result.Add((i, j, l));
            }
        return result;
    }

    private static IEnumerable<ISuffixTree> Trees(string text)
    {
        yield return global::SuffixTree.SuffixTree.Build(text);
        yield return PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text));
    }

    [Test]
    public void MaximalPairs_EqualBruteForce_MinLengthNesting_Relabelling()
    {
        var rng = new Random(260_001);
        Func<char, bool>?[] predicates = { null, c => c == 'N', c => "ACGT".IndexOf(c) < 0 };
        for (int trial = 0; trial < 120; trial++)
        {
            string alphabet = (trial % 3) switch { 0 => "AC", 1 => "ACGT", _ => "ACGTN" };
            string text = RandomString(rng, rng.Next(0, 70), alphabet);
            int minLength = rng.Next(1, 5);
            var unique = predicates[rng.Next(predicates.Length)];
            var expected = BruteForcePairs(text, minLength, unique);
            string msg = $"trial {trial} {text} min={minLength}";
            foreach (var tree in Trees(text))
            {
                try
                {
                    var pairs = tree.FindMaximalRepeatedPairs(minLength, unique);
                    Assert.That(pairs.Select(p => (p.FirstPosition, p.SecondPosition, p.Length)), Is.EqualTo(expected), msg);
                    int d = rng.Next(1, 4);
                    Assert.That(tree.FindMaximalRepeatedPairs(minLength + d, unique),
                        Is.EqualTo(pairs.Where(p => p.Length >= minLength + d).ToList()), msg + " nesting");
                }
                finally
                {
                    (tree as IDisposable)?.Dispose();
                }
            }

            if (unique is null)
            {
                // Bijection A→T, C→G, G→C, T→A, N→N (complement without reversal).
                string relabelled = new(text.Select(c => c switch { 'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', _ => c }).ToArray());
                Assert.That(global::SuffixTree.SuffixTree.Build(relabelled).FindMaximalRepeatedPairs(minLength),
                    Is.EqualTo(global::SuffixTree.SuffixTree.Build(text).FindMaximalRepeatedPairs(minLength)), msg + " relabel");
            }
        }
    }

    [Test]
    public void MaximalPairs_UniquePredicate_EqualsFreshSymbolSubstitution()
    {
        var rng = new Random(260_002);
        for (int trial = 0; trial < 100; trial++)
        {
            string text = RandomString(rng, rng.Next(0, 120), "ACGTN");
            int minLength = rng.Next(1, 5);
            int fresh = 0xE000; // private-use code points, one per N occurrence
            string substituted = new(text.Select(c => c == 'N' ? (char)fresh++ : c).ToArray());
            var expected = global::SuffixTree.SuffixTree.Build(substituted).FindMaximalRepeatedPairs(minLength);
            using var persistent = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text));
            Assert.That(global::SuffixTree.SuffixTree.Build(text).FindMaximalRepeatedPairs(minLength, c => c == 'N'), Is.EqualTo(expected), text);
            Assert.That(persistent.FindMaximalRepeatedPairs(minLength, c => c == 'N'), Is.EqualTo(expected), text);
        }
    }

    private static int BruteLcsLength(string a, string b)
    {
        int best = 0;
        var prev = new int[b.Length + 1];
        for (int i = 1; i <= a.Length; i++)
        {
            var cur = new int[b.Length + 1];
            for (int j = 1; j <= b.Length; j++)
                if (a[i - 1] == b[j - 1])
                {
                    cur[j] = prev[j - 1] + 1;
                    best = Math.Max(best, cur[j]);
                }
            prev = cur;
        }
        return best;
    }

    private static List<int> Occurrences(string s, string w)
    {
        var r = new List<int>();
        for (int i = 0; i + w.Length <= s.Length; i++)
            if (string.CompareOrdinal(s, i, w, 0, w.Length) == 0) r.Add(i);
        return r;
    }

    [Test]
    public void AllDistinctLcs_SwapSymmetric_LengthAndPositionsByBruteForce()
    {
        var rng = new Random(260_003);
        for (int trial = 0; trial < 150; trial++)
        {
            string alphabet = trial % 2 == 0 ? "ACG" : "ACGT";
            string a = RandomString(rng, rng.Next(0, 60), alphabet);
            string b = RandomString(rng, rng.Next(0, 60), alphabet);
            int lcs = BruteLcsLength(a, b);
            string msg = $"trial {trial} {a} | {b}";

            using var pa = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(a));
            var ab = pa.FindAllDistinctLongestCommonSubstrings(b);
            var ba = global::SuffixTree.SuffixTree.Build(b).FindAllDistinctLongestCommonSubstrings(a);
            Assert.That(ab.Select(x => x.Substring), Is.EquivalentTo(ba.Select(x => x.Substring)), msg);
            Assert.That(ab.Select(x => x.Substring), Is.Unique, msg);
            if (lcs == 0)
                Assert.That(ab, Is.Empty, msg);
            foreach (var x in ab)
            {
                Assert.That(x.Substring.Length, Is.EqualTo(lcs), msg);
                Assert.That(x.PositionsInText, Is.EqualTo(Occurrences(a, x.Substring)), msg);
                Assert.That(x.PositionsInOther, Is.EqualTo(Occurrences(b, x.Substring)), msg);
            }
            // Every common substring of length lcs is reported.
            var all = Enumerable.Range(0, Math.Max(0, a.Length - lcs + 1)).Select(i => a.Substring(i, lcs))
                .Where(w => lcs > 0 && b.Contains(w, StringComparison.Ordinal)).Distinct();
            Assert.That(ab.Select(x => x.Substring), Is.EquivalentTo(all), msg);
        }
    }

    [Test]
    public void KCommonSubstrings_Monotone_PermutationInvariant_BruteForce()
    {
        var rng = new Random(260_004);
        Func<string, ISuffixTree> persistent = s => PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(s));
        for (int trial = 0; trial < 60; trial++)
        {
            int k = rng.Next(1, 6);
            var texts = Enumerable.Range(0, k).Select(_ => RandomString(rng, rng.Next(0, 40), trial % 2 == 0 ? "AC" : "ACGT")).ToArray();
            string msg = $"trial {trial} [{string.Join(",", texts)}]";
            int[] l = global::SuffixTree.SuffixTree.LongestCommonSubstringLengthsBySupport(texts);
            Assert.That(l[1], Is.EqualTo(texts.Max(t => t.Length)), msg);
            for (int q = 2; q <= k; q++)
                Assert.That(l[q], Is.LessThanOrEqualTo(l[q - 1]), msg + $" q={q}");

            var perm = texts.OrderBy(_ => rng.Next()).ToArray();
            Assert.That(SuffixTreeAlgorithms.LongestCommonSubstringLengthsBySupport(perm, persistent), Is.EqualTo(l), msg + " perm");

            // Brute force: support of every distinct substring.
            var support = new Dictionary<string, int>(StringComparer.Ordinal);
            foreach (var t in texts)
            {
                var seen = new HashSet<string>(StringComparer.Ordinal);
                for (int i = 0; i < t.Length; i++)
                    for (int j = i + 1; j <= t.Length; j++)
                        seen.Add(t.Substring(i, j - i));
                foreach (var w in seen)
                    support[w] = support.GetValueOrDefault(w) + 1;
            }
            for (int q = 1; q <= k; q++)
            {
                var result = global::SuffixTree.SuffixTree.FindLongestCommonSubstrings(texts, q);
                int bruteL = support.Where(p => p.Value >= q).Select(p => p.Key.Length).DefaultIfEmpty(0).Max();
                Assert.That(l[q], Is.EqualTo(bruteL), msg + $" l({q})");
                var expected = support.Where(p => p.Value >= q && p.Key.Length == bruteL && bruteL > 0).Select(p => p.Key);
                Assert.That(result, Is.EquivalentTo(expected), msg + $" q={q}");
                Assert.That(SuffixTreeAlgorithms.FindLongestCommonSubstrings(perm, q, persistent), Is.EqualTo(result), msg + $" perm q={q}");
            }
        }
    }
}
