using SuffixTree.Esa.Tests.TestSupport;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Tree;

/// <summary>The in-memory ESA (work package 1b.1): arrays, header statistics, n = 0.</summary>
[TestFixture]
public class EsaIndexTests
{
    private static IEnumerable<string> Corpus()
    {
        foreach (var (sigma, maxLen) in new[] { (2, 8), (3, 5), (4, 4) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
                yield return s;
        var rng = new SplitMix64(0x1B04);
        for (int i = 0; i < 12; i++)
        {
            int n = rng.NextInt(300);
            ulong seed = rng.Next();
            foreach (var (_, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 5 == 0))
                yield return text;
        }
    }

    [Test]
    public void Arrays_EqualTheNaiveOracles()
    {
        foreach (string s in Corpus())
        {
            EsaIndex index = EsaIndex.Build(s);
            int[] sa = NaiveSuffixArray.Build(s);
            Assert.That(index.SuffixArray, Is.EqualTo(sa), $"'{s}'");
            Assert.That(index.Lcp.ToArray(), Is.EqualTo(NaiveLcp.Build(s, sa)), $"'{s}'");
            Assert.That(index.Text, Is.SameAs(s));
            Assert.That(index.Length, Is.EqualTo(s.Length));
        }
    }

    [Test]
    public void Statistics_EqualTheUkkonenTree()
    {
        foreach (string s in Corpus())
        {
            UkkonenTree tree = UkkonenTree.Build(s);
            EsaStatistics st = EsaIndex.Build(s).Statistics;
            Assert.That(st.NodeCount, Is.EqualTo(tree.NodeCount), $"NodeCount '{s}'");
            Assert.That(st.LeafCount, Is.EqualTo(tree.LeafCount), $"LeafCount '{s}'");
            Assert.That(st.MaxDepth, Is.EqualTo(tree.MaxDepth), $"MaxDepth '{s}'");
            Assert.That(st.Length, Is.EqualTo(s.Length), s);
        }
    }

    [Test]
    public void Statistics_NodeCountIsRootPlusInternalPlusLeaves()
    {
        EsaStatistics st = EsaIndex.Build("mississippi").Statistics;
        // mississippi: 11 suffixes + the empty suffix = 12 leaves; internal nodes besides the root: i, issi, ssi, si, p, ... counted by the scan
        Assert.That(st.NodeCount, Is.EqualTo(1 + st.InternalCount + 12));
        Assert.That(st.LeafCount, Is.EqualTo(11));
        Assert.That(st.MaxDepth, Is.EqualTo(11));
    }

    [Test]
    public void EmptyText_HasOnlyTheRoot()
    {
        EsaIndex index = EsaIndex.Build("");
        Assert.That(index.Length, Is.EqualTo(0));
        Assert.That(index.SuffixArray, Is.EqualTo(new[] { 0 }));
        Assert.That(index.Statistics, Is.EqualTo(new EsaStatistics(0, 1, 0, 0, 0, 0, 0, 0)));
    }

    [Test]
    public void LongestRepeatedSubstring_DepthAndRankDescribeADeepestInternalNode()
    {
        foreach (string s in Corpus().Where(s => s.Length > 0))
        {
            EsaIndex index = EsaIndex.Build(s);
            EsaStatistics st = index.Statistics;
            UkkonenTree tree = UkkonenTree.Build(s);
            string lrs = tree.LongestRepeatedSubstring();
            Assert.That(st.LrsDepth, Is.EqualTo(lrs.Length), $"'{s}'");
            if (st.LrsDepth > 0)
            {
                // the rank is the left boundary of an lcp-interval of that depth: the suffix at that rank starts with the repeat
                int count = BruteForce.Occurrences(s, s.Substring(index.SuffixArray[st.LrsRank], st.LrsDepth)).Count;
                Assert.That(count, Is.GreaterThanOrEqualTo(2), $"'{s}'");
            }
        }
    }

    [Test]
    public void Build_NullText_Throws()
    {
        Assert.That(() => EsaIndex.Build(null!), Throws.TypeOf<ArgumentNullException>());
    }
}
