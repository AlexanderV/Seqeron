using SuffixTree.Esa.Tests.TestSupport;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Oracles;

/// <summary>1a.1: the test-only definitional oracles (naive SA/LCP, Kasai) and their cross-checks.</summary>
[TestFixture]
public class NaiveOracleTests
{
    [Test]
    public void NaiveSuffixArray_Banana()
    {
        Assert.That(NaiveSuffixArray.Build("banana"), Is.EqualTo(new[] { 6, 5, 3, 1, 0, 4, 2 }));
    }

    [Test]
    public void NaiveSuffixArray_TinyTexts()
    {
        Assert.That(NaiveSuffixArray.Build(""), Is.EqualTo(new[] { 0 }));
        Assert.That(NaiveSuffixArray.Build("x"), Is.EqualTo(new[] { 1, 0 }));
        Assert.That(NaiveSuffixArray.Build("ba"), Is.EqualTo(new[] { 2, 1, 0 }));
        Assert.That(NaiveSuffixArray.Build("ab"), Is.EqualTo(new[] { 2, 0, 1 }));
        Assert.That(NaiveSuffixArray.Build("aaaa"), Is.EqualTo(new[] { 4, 3, 2, 1, 0 }));
    }

    [Test]
    public void NaiveSuffixArray_UsesOrdinalCodeUnitOrder()
    {
        // U+0000 < 'a' < U+D800 (lone surrogate) < U+FFFF in code-unit order, regardless of culture.
        string s = "￿a\uD800\0";
        Assert.That(NaiveSuffixArray.Build(s), Is.EqualTo(new[] { 4, 3, 1, 2, 0 }));
        // A proper prefix sorts before its extensions (the virtual terminator is the smallest symbol).
        Assert.That(NaiveSuffixArray.Build("abab"), Is.EqualTo(new[] { 4, 2, 0, 3, 1 }));
    }

    [Test]
    public void NaiveSuffixArray_IsAStrictlyIncreasingPermutation()
    {
        foreach (var (name, text) in StringFamilies.Adversarial(120, 5))
        {
            int[] sa = NaiveSuffixArray.Build(text);
            Assert.That(sa, Has.Length.EqualTo(text.Length + 1), name);
            Assert.That(sa[0], Is.EqualTo(text.Length), name);
            Assert.That(sa.OrderBy(x => x), Is.EqualTo(Enumerable.Range(0, text.Length + 1)), name);
            for (int r = 1; r < sa.Length; r++)
                Assert.That(text.AsSpan(sa[r - 1]).SequenceCompareTo(text.AsSpan(sa[r])), Is.LessThan(0), $"{name} rank {r}");
        }
    }

    [Test]
    public void NaiveSuffixArray_EqualsUkkonenSuffixOrder_On10000RandomAndAdversarialTexts()
    {
        var rng = new SplitMix64(20261103);
        int[] sigmas = { 1, 2, 3, 4, 20, 256 };
        for (int i = 0; i < 10_000; i++)
        {
            int n = 1 + rng.NextInt(500);
            ulong seed = rng.Next();
            string text;
            if ((i & 1) == 0)
            {
                text = StringFamilies.Random(n, sigmas[rng.NextInt(sigmas.Length)], seed);
            }
            else
            {
                var all = StringFamilies.Adversarial(n, seed).ToList();
                text = all[rng.NextInt(all.Count)].Text;
            }

            int[] sa = NaiveSuffixArray.Build(text);
            var expected = UkkonenTree.Build(text).EnumerateSuffixes().ToList();
            Assert.That(expected, Has.Count.EqualTo(text.Length), $"case {i} n={n}");
            for (int r = 1; r <= text.Length; r++)
            {
                if (!string.Equals(text.Substring(sa[r]), expected[r - 1], StringComparison.Ordinal))
                    Assert.Fail($"case {i} seed={seed} n={n} rank {r}: naive suffix at {sa[r]} differs from Ukkonen order");
            }
        }
    }

    [Test]
    public void NaiveLcp_Banana()
    {
        int[] sa = NaiveSuffixArray.Build("banana");
        Assert.That(NaiveLcp.Build("banana", sa), Is.EqualTo(new[] { 0, 0, 1, 3, 0, 0, 2 }));
    }

    [Test]
    public void NaiveLcp_TinyTexts()
    {
        Assert.That(NaiveLcp.Build("", NaiveSuffixArray.Build("")), Is.EqualTo(new[] { 0 }));
        Assert.That(NaiveLcp.Build("x", NaiveSuffixArray.Build("x")), Is.EqualTo(new[] { 0, 0 }));
        // a^n: LCP[r] = r - 1
        Assert.That(NaiveLcp.Build("aaaa", NaiveSuffixArray.Build("aaaa")), Is.EqualTo(new[] { 0, 0, 1, 2, 3 }));
    }

    [Test]
    public void Kasai_Banana_AndTinyTexts()
    {
        Assert.That(KasaiLcp.Build("banana", NaiveSuffixArray.Build("banana")), Is.EqualTo(new[] { 0, 0, 1, 3, 0, 0, 2 }));
        Assert.That(KasaiLcp.Build("", NaiveSuffixArray.Build("")), Is.EqualTo(new[] { 0 }));
        Assert.That(KasaiLcp.Build("x", NaiveSuffixArray.Build("x")), Is.EqualTo(new[] { 0, 0 }));
    }

    [Test]
    public void Kasai_EqualsNaiveLcp_OnAllSmallStrings()
    {
        foreach (var (sigma, maxLen) in new[] { (2, 10), (3, 6), (4, 5) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
            {
                int[] sa = NaiveSuffixArray.Build(s);
                if (!NaiveLcp.Build(s, sa).AsSpan().SequenceEqual(KasaiLcp.Build(s, sa)))
                    Assert.Fail($"Kasai differs from naive LCP on '{s}'");
            }
    }

    [Test]
    public void Kasai_EqualsNaiveLcp_OnRandomAndAdversarialTexts()
    {
        var rng = new SplitMix64(77);
        for (int i = 0; i < 400; i++)
        {
            int n = rng.NextInt(600);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed))
            {
                int[] sa = NaiveSuffixArray.Build(text);
                if (!NaiveLcp.Build(text, sa).AsSpan().SequenceEqual(KasaiLcp.Build(text, sa)))
                    Assert.Fail($"Kasai differs from naive LCP: family {name} seed {seed} n {n}");
            }
        }
    }

    [Test]
    public void ExhaustiveStrings_EnumeratesEveryStringOnce()
    {
        var all = ExhaustiveStrings.All(2, 3).ToList();
        Assert.That(all, Is.EqualTo(new[] { "", "a", "b", "aa", "ab", "ba", "bb", "aaa", "aab", "aba", "abb", "baa", "bab", "bba", "bbb" }));
        Assert.That(ExhaustiveStrings.All(3, 4).Count(), Is.EqualTo(1 + 3 + 9 + 27 + 81));
        Assert.That(ExhaustiveStrings.All(3, 4).Distinct().Count(), Is.EqualTo(121));
    }
}
