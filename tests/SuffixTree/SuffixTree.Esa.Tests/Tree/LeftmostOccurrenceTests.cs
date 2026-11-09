using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Tree;

/// <summary>
/// Hypothesis H-OCC (ADR-001 section 4) on the real engine: the edge of every internal node records the leftmost occurrence of its string,
/// <c>min SA[lb..rb]</c>. <see cref="LeftmostOccurrence"/> computes it for all internal intervals in one bottom-up scan; this fixture checks it against
/// the definition (the minimum over the suffix-array range). That the Ukkonen and v6 trees agree with the rule is checked by the Traverse streams in
/// <see cref="EsaTraverseTests"/>.
/// </summary>
[TestFixture]
public class LeftmostOccurrenceTests
{
    private static void AssertMatchesDefinition(string text, string context)
    {
        EsaIndex index = EsaIndex.Build(text);
        int[] occ = LeftmostOccurrence.Compute(index.SuffixArray, index.Lcp);
        foreach (BruteInterval x in BruteLcpTree.Intervals(index.Lcp.ToArray()))
        {
            int expected = int.MaxValue;
            for (int r = x.Lb; r <= x.Rb; r++) expected = Math.Min(expected, index.SuffixArray[r]);
            Assert.That(occ[x.LIndices[0]], Is.EqualTo(expected), $"{context} interval [{x.Lb}..{x.Rb}]");
        }
    }

    [Test]
    public void Banana_EveryInternalNode()
    {
        AssertMatchesDefinition("banana", "banana");
        AssertMatchesDefinition("mississippi", "mississippi");
        AssertMatchesDefinition("abracadabra", "abracadabra");
    }

    [Test]
    public void AllSmallStrings()
    {
        foreach (var (sigma, maxLen) in new[] { (2, 10), (3, 6), (4, 5) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen).Where(s => s.Length > 0))
                AssertMatchesDefinition(s, $"'{s}'");
    }

    [Test]
    public void AdversarialAndRandomTexts()
    {
        var rng = new SplitMix64(0x1B05);
        for (int i = 0; i < 30; i++)
        {
            int n = 1 + rng.NextInt(300);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 3 == 0))
                AssertMatchesDefinition(text, $"{name} n {n}");
        }
    }

    [Test]
    public void EmptyText_HasNoInternalNode()
    {
        EsaIndex index = EsaIndex.Build("");
        Assert.That(LeftmostOccurrence.Compute(index.SuffixArray, index.Lcp), Has.Length.EqualTo(1));
    }
}
