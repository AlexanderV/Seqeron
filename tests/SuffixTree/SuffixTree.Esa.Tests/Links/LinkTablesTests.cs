using SuffixTree.Esa.Bench.Links;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Links;

/// <summary>ISA and Psi (work package 1b.4a): the definitions, the monotonicity that identity F1 of ADR-001 section 5 rests on, and the node enumeration used by the measurements.</summary>
[TestFixture]
public class LinkTablesTests
{
    private static void CheckText(string text)
    {
        int n = text.Length;
        int[] sa = NaiveSuffixArray.Build(text);
        int[] isa = LinkTables.BuildIsa(sa);
        int[] psi = LinkTables.BuildPsi(sa, isa);
        if (isa.Length != n + 1 || psi.Length != n + 1) Assert.Fail($"'{text}': lengths {isa.Length} / {psi.Length}");
        for (int r = 0; r <= n; r++)
        {
            if (isa[sa[r]] != r) Assert.Fail($"'{text}': isa[sa[{r}]] = {isa[sa[r]]}");
            if (sa[isa[r]] != r) Assert.Fail($"'{text}': sa[isa[{r}]] = {sa[isa[r]]}");
        }
        if (isa[n] != 0) Assert.Fail($"'{text}': the empty suffix must have rank 0");
        if (psi[0] != -1) Assert.Fail($"'{text}': psi[0] = {psi[0]}");
        for (int r = 1; r <= n; r++)
        {
            // Psi(r) is the rank of the suffix that starts one symbol later (the empty suffix, rank 0, after the last symbol)
            int expected = isa[sa[r] + 1];
            if (psi[r] != expected) Assert.Fail($"'{text}': psi[{r}] = {psi[r]}, expected {expected}");
            if (sa[psi[r]] != sa[r] + 1) Assert.Fail($"'{text}': psi[{r}] does not point at the successor suffix");
        }
        // monotone on every run of suffixes with a common first symbol (the order of the suffixes is the order of their tails)
        for (int r = 2; r <= n; r++)
            if (text[sa[r - 1]] == text[sa[r]] && psi[r - 1] >= psi[r]) Assert.Fail($"'{text}': psi is not increasing between ranks {r - 1} and {r}");
    }

    [Test]
    public void Definitions_OnAllSmallStrings()
    {
        foreach (string s in ExhaustiveStrings.All(2, 10)) CheckText(s);
        foreach (string s in ExhaustiveStrings.All(3, 6)) CheckText(s);
        foreach (string s in ExhaustiveStrings.All(4, 4)) CheckText(s);
    }

    [Test]
    public void Definitions_OnAdversarialFamilies()
    {
        foreach (int n in new[] { 1, 2, 3, 17, 100, 301 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 7))
            {
                try
                {
                    CheckText(text);
                }
                catch (AssertionException e)
                {
                    Assert.Fail($"{name} n {n}: {e.Message}");
                }
            }
    }

    [Test]
    public void ArraysAreFourBytesPerEntry_OnTheRealIndex()
    {
        EsaIndex esa = EsaIndex.Build(StringFamilies.RandomDna(5000, 3));
        int[] isa = LinkTables.BuildIsa(esa.SuffixArray);
        int[] psi = LinkTables.BuildPsi(esa.SuffixArray, isa);
        Assert.That(isa, Has.Length.EqualTo(5001));
        Assert.That(psi, Has.Length.EqualTo(5001));
        for (int r = 1; r <= 5000; r++) Assert.That(psi[r], Is.EqualTo(isa[esa.SuffixArray[r] + 1]));
    }
}

/// <summary>The enumeration of internal nodes with their depths, the deterministic sample and the largest-key selection.</summary>
[TestFixture]
public class LinkNodeSetTests
{
    [Test]
    public void Enumerate_EqualsTheDefinitionalIntervalTree()
    {
        foreach (string text in ExhaustiveStrings.All(2, 9).Concat(ExhaustiveStrings.All(3, 5)).Concat(StringFamilies.Adversarial(200, 5).Select(f => f.Text)))
        {
            int[] sa = NaiveSuffixArray.Build(text);
            int[] lcp = NaiveLcp.Build(text, sa);
            var brute = BruteLcpTree.Intervals(lcp).Where(iv => iv.Ell >= 1).Select(iv => (iv.Lb, iv.Rb, iv.Ell)).OrderBy(x => x.Lb).ThenByDescending(x => x.Rb).ToList();
            var set = LinkNodeSet.Enumerate(CompactLcp.FromDense(lcp));
            var actual = Enumerable.Range(0, set.Count).Select(i => (set.Lb[i], set.Rb[i], set.Depth[i])).OrderBy(x => x.Item1).ThenByDescending(x => x.Item2).ToList();
            if (!brute.SequenceEqual(actual)) Assert.Fail($"'{(text.Length > 40 ? text[..40] + "..." : text)}': internal nodes differ ({brute.Count} against {actual.Count})");

            // node depth: the number of internal ancestors including the node itself (a child of the root has 1)
            for (int i = 0; i < set.Count; i++)
            {
                int ancestors = actual.Count(x => x.Item1 <= set.Lb[i] && set.Rb[i] <= x.Item2);
                if (set.NodeDepth[i] != ancestors) Assert.Fail($"'{text}': node ({set.Lb[i]}, {set.Rb[i]}) has NodeDepth {set.NodeDepth[i]}, expected {ancestors}");
            }
        }
    }

    [Test]
    public void Enumerate_CountMatchesTheIndexStatistics()
    {
        foreach (string text in new[] { StringFamilies.RandomDna(20000, 1), StringFamilies.Fibonacci(20000), StringFamilies.AllSame(3000), StringFamilies.RandomProtein(9000, 2) })
        {
            EsaIndex esa = EsaIndex.Build(text);
            Assert.That(LinkNodeSet.Enumerate(esa.Lcp).Count, Is.EqualTo(esa.Statistics.InternalCount));
        }
    }

    [Test]
    public void Node_CarriesTheFirstSuffixOfTheInterval()
    {
        EsaIndex esa = EsaIndex.Build(StringFamilies.RandomDna(500, 4));
        var set = LinkNodeSet.Enumerate(esa.Lcp);
        for (int i = 0; i < set.Count; i++)
        {
            LinkNode v = set.Node(i, esa.SuffixArray);
            Assert.That((v.Lb, v.Rb, v.Depth, v.Pos), Is.EqualTo((set.Lb[i], set.Rb[i], set.Depth[i], esa.SuffixArray[set.Lb[i]])));
        }
    }

    [Test]
    public void Sample_IsDeterministicShuffledAndWithoutRepetition()
    {
        EsaIndex esa = EsaIndex.Build(StringFamilies.RandomDna(20000, 6));
        var set = LinkNodeSet.Enumerate(esa.Lcp);
        int[] a = set.Sample(1000, 42), b = set.Sample(1000, 42), c = set.Sample(1000, 43);
        Assert.That(a, Is.EqualTo(b));
        Assert.That(a, Is.Not.EqualTo(c));
        Assert.That(a.Distinct().Count(), Is.EqualTo(1000));
        Assert.That(a.All(i => (uint)i < (uint)set.Count));
        Assert.That(a, Is.Not.Ordered, "the sample is shuffled, not in suffix-array order");
        int[] all = set.Sample(set.Count + 50, 1);
        Assert.That(all.OrderBy(x => x), Is.EqualTo(Enumerable.Range(0, set.Count)), "fewer nodes than asked: all of them");
    }

    [Test]
    public void Largest_ReturnsTheKLargestKeysLargestFirst()
    {
        var rng = new SplitMix64(9);
        for (int rep = 0; rep < 50; rep++)
        {
            int[] keys = Enumerable.Range(0, 1 + rng.NextInt(300)).Select(_ => rng.NextInt(40)).ToArray();
            int k = 1 + rng.NextInt(20);
            int[] got = LinkNodeSet.Largest(keys, k);
            int[] expected = Enumerable.Range(0, keys.Length).OrderByDescending(i => keys[i]).ThenBy(i => i).Take(k).ToArray();
            Assert.That(got, Is.EqualTo(expected), $"rep {rep}");
        }
    }
}
