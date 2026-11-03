using SuffixTree.Esa.Tests.TestSupport;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Oracles;

/// <summary>
/// Hypothesis H-OCC (ADR-001 section 4; review MIN-3), tested test-only right after the naive oracles: a naive lcp-interval
/// tree built from the naive SA and LCP, with the edge of every internal node labelled by the LEFTMOST occurrence
/// (min SA over the interval), produces exactly the <c>Traverse</c> event stream of the Ukkonen tree. This is the content
/// of <c>SuffixTreeSerializer.CalculateLogicalHash</c> (text bytes + the event stream).
/// </summary>
[TestFixture]
public class HoccModelTests
{
    private static void AssertSameStream(string text, string context)
    {
        var expected = new TraverseRecorder();
        UkkonenTree.Build(text).Traverse(expected);
        var actual = new TraverseRecorder();
        NaiveEsaModel.EmitTraverse(text, actual);
        if (!expected.Events.SequenceEqual(actual.Events))
        {
            int i = 0;
            while (i < expected.Events.Count && i < actual.Events.Count && expected.Events[i] == actual.Events[i]) i++;
            Assert.Fail($"{context}: first difference at event {i}: Ukkonen '{(i < expected.Events.Count ? expected.Events[i] : "<end>")}' vs model '{(i < actual.Events.Count ? actual.Events[i] : "<end>")}' (text length {text.Length})");
        }
    }

    [Test]
    public void Banana_StreamMatches()
    {
        AssertSameStream("banana", "banana");
        AssertSameStream("mississippi", "mississippi");
        AssertSameStream("abracadabra", "abracadabra");
    }

    [Test]
    public void TinyTexts_StreamMatches()
    {
        foreach (string s in new[] { "", "a", "aa", "ab", "ba", "aaa", "aba" })
            AssertSameStream(s, $"'{s}'");
    }

    [Test]
    public void AllSmallStrings_StreamMatches()
    {
        foreach (var (sigma, maxLen) in new[] { (2, 11), (3, 7), (4, 5) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
                AssertSameStream(s, $"'{s}'");
    }

    [Test]
    public void RandomAndAdversarialTexts_StreamMatches()
    {
        var rng = new SplitMix64(31415);
        for (int i = 0; i < 150; i++)
        {
            int n = rng.NextInt(300);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed))
                AssertSameStream(text, $"family {name} seed {seed} n {n}");
        }
    }
}
