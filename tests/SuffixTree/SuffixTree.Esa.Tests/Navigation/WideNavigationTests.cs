using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.Search;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Navigation;

/// <summary>
/// Work package S2.5, the navigation side at 64 bits (written before the engine: against the engine of S2.4 the cases over 5-byte entries are red, they end in the named refusal of the 32-bit pack). Every operation that needs the
/// navigation pack (<c>TrySuffixLink</c>, <c>TryParent</c>, <c>Lcp</c>, <c>Lca</c>, <c>RankOf</c>, <c>TryExtendLeft</c>) is applied to every explicit node of the tree and to the ranks and positions of the text, over every layout of the
/// matrix of <see cref="WideImage"/> (4- and 5-byte entries, the LCP array of 8, 16, 32 and 64 bits, packed and tagged CHILD) in both access modes, and compared with the same operation of the index that was built in memory.
/// </summary>
[TestFixture]
public class WideNavigationTests
{
    private static readonly EsaAccessMode[] Modes = [EsaAccessMode.Checked, EsaAccessMode.Trusted];

    public static IEnumerable<string> Texts()
    {
        var rng = new SplitMix64(0x25C0DE);
        yield return "mississippi";
        yield return "abracadabra";
        yield return "aaaaaaaaaaaaaaaaaaaa";
        yield return "abababababababababab";
        yield return "banana";
        yield return "a";
        yield return RandomText(rng, 300, 4);
        yield return RandomText(rng, 200, 2);
        yield return RandomText(rng, 150, 20);
        string block = RandomText(rng, 40, 4);
        yield return RandomText(rng, 60, 4) + block + RandomText(rng, 60, 4) + block + new string('A', 200);
    }

    private static string RandomText(SplitMix64 rng, int n, int sigma)
    {
        var chars = new char[n];
        for (int i = 0; i < n; i++) chars[i] = (char)('A' + rng.NextInt(sigma));
        return new string(chars);
    }

    private static void ForEachLayout(string text, Action<string, EsaSuffixIndex, EsaSuffixIndex> check)
    {
        var memory = new EsaSuffixIndex(EsaIndex.Build(text));
        LcpChildCase c = LcpChildCase.Of(text);
        char[] alphabet = [.. text.Distinct().Order()];
        foreach ((int width, LcpKind kind, bool tagged) in WideImage.Layouts)
        {
            using var image = new WideImage(c, width, kind, tagged);
            foreach (EsaAccessMode mode in Modes)
            {
                var wide = new EsaSuffixIndex(image.Index(mode == EsaAccessMode.Checked, alphabet));
                check($"width {width} {kind}{(tagged ? " tagged" : string.Empty)} {mode} '{(text.Length <= 16 ? text : text[..16] + "...")}'", memory, wide);
            }
        }
    }

    /// <summary>Every explicit node of the tree (the root included), breadth first, read through the contract.</summary>
    internal static List<IndexLocus> Nodes(EsaSuffixIndex idx)
    {
        var all = new List<IndexLocus> { idx.Root };
        var buffer = new IndexLocus[64];
        for (int i = 0; i < all.Count; i++)
        {
            int count = idx.GetChildren(all[i], buffer);
            Assert.That(count, Is.LessThanOrEqualTo(buffer.Length));
            for (int k = 0; k < count; k++) all.Add(buffer[k]);
        }
        return all;
    }

    /// <summary>Every navigation operation of <paramref name="w"/> against <paramref name="m"/> (the same text): the links, the parents and the left extensions of every node, the arrays and range queries on a grid of ranks.</summary>
    internal static void CompareAll(string label, EsaSuffixIndex m, EsaSuffixIndex w, string text)
    {
        List<IndexLocus> nodes = Nodes(m);
        Assert.That(Nodes(w), Is.EqualTo(nodes), $"{label} the nodes");
        int[] symbols = [.. text.Distinct().Select(ch => (int)ch), 'z'];
        foreach (IndexLocus node in nodes)
        {
            Assert.That(w.TrySuffixLink(node, out IndexLocus link), Is.EqualTo(m.TrySuffixLink(node, out IndexLocus expectedLink)), $"{label} TrySuffixLink {node}");
            Assert.That(link, Is.EqualTo(expectedLink), $"{label} link of {node}");
            Assert.That(w.TryParent(node, out IndexLocus parent), Is.EqualTo(m.TryParent(node, out IndexLocus expectedParent)), $"{label} TryParent {node}");
            Assert.That(parent, Is.EqualTo(expectedParent), $"{label} parent of {node}");
            foreach (int symbol in symbols)
            {
                Assert.That(w.TryExtendLeft(node, symbol, out IndexLocus extended), Is.EqualTo(m.TryExtendLeft(node, symbol, out IndexLocus expectedExtended)), $"{label} TryExtendLeft {node} {symbol}");
                Assert.That(extended, Is.EqualTo(expectedExtended), $"{label} extend left {node} {symbol}");
            }
        }
        int step = Math.Max(1, text.Length / 40);
        for (long a = 0; a <= text.Length; a += step)
        {
            Assert.That(w.RankOf(a), Is.EqualTo(m.RankOf(a)), $"{label} ISA[{a}]");
            for (long b = 0; b <= text.Length; b += step)
            {
                Assert.That(w.Lcp(a, b), Is.EqualTo(m.Lcp(a, b)), $"{label} Lcp({a}, {b})");
                Assert.That(w.Lca(a, b), Is.EqualTo(m.Lca(a, b)), $"{label} Lca({a}, {b})");
            }
        }
    }

    [Test]
    public void TheSuffixLink_OfEveryNode_EqualsTheMemoryIndex_AtEveryLayout([ValueSource(nameof(Texts))] string text)
    {
        ForEachLayout(text, (label, m, w) =>
        {
            List<IndexLocus> nodes = Nodes(m);
            Assert.That(Nodes(w), Is.EqualTo(nodes), $"{label} the nodes");
            foreach (IndexLocus node in nodes)
            {
                Assert.That(w.TrySuffixLink(node, out IndexLocus link), Is.EqualTo(m.TrySuffixLink(node, out IndexLocus expected)), $"{label} TrySuffixLink {node}");
                Assert.That(link, Is.EqualTo(expected), $"{label} link of {node}");
            }
        });
    }

    [Test]
    public void TheParent_OfEveryNodeAndOfTheFirstSymbolPoint_EqualsTheMemoryIndex_AtEveryLayout([ValueSource(nameof(Texts))] string text)
    {
        ForEachLayout(text, (label, m, w) =>
        {
            foreach (IndexLocus node in Nodes(m))
            {
                Assert.That(w.TryParent(node, out IndexLocus parent), Is.EqualTo(m.TryParent(node, out IndexLocus expected)), $"{label} TryParent {node}");
                Assert.That(parent, Is.EqualTo(expected), $"{label} parent of {node}");
            }
        });
    }

    [Test]
    public void TheLcpAndTheLca_OfEveryPairOfRanks_AndTheRankOfEveryPosition_EqualTheMemoryIndex_AtEveryLayout([ValueSource(nameof(Texts))] string text)
    {
        ForEachLayout(text, (label, m, w) =>
        {
            int step = Math.Max(1, text.Length / 40);
            for (long a = 0; a <= text.Length; a += step)
            {
                Assert.That(w.RankOf(a), Is.EqualTo(m.RankOf(a)), $"{label} ISA[{a}]");
                for (long b = 0; b <= text.Length; b += step)
                {
                    Assert.That(w.Lcp(a, b), Is.EqualTo(m.Lcp(a, b)), $"{label} Lcp({a}, {b})");
                    Assert.That(w.Lca(a, b), Is.EqualTo(m.Lca(a, b)), $"{label} Lca({a}, {b})");
                }
            }
        });
    }

    [Test]
    public void TheLeftExtension_OfEveryNode_ByEverySymbolOfTheAlphabet_EqualsTheMemoryIndex_AtEveryLayout([ValueSource(nameof(Texts))] string text)
    {
        int[] symbols = [.. text.Distinct().Select(ch => (int)ch), 'z'];
        ForEachLayout(text, (label, m, w) =>
        {
            foreach (IndexLocus node in Nodes(m))
            {
                foreach (int symbol in symbols)
                {
                    Assert.That(w.TryExtendLeft(node, symbol, out IndexLocus extended), Is.EqualTo(m.TryExtendLeft(node, symbol, out IndexLocus expected)), $"{label} TryExtendLeft {node} {symbol}");
                    Assert.That(extended, Is.EqualTo(expected), $"{label} extend left {node} {symbol}");
                }
            }
        });
    }
}
