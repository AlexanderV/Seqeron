using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Tree;

/// <summary>
/// The three child-table layouts of prototype P3 (ADR-001 section 5, work package 1b.1): the plain three-array reference, the single-array packing
/// whose field kinds are told apart by LCP comparisons (Abouelhoda, Kurtz, Ohlebusch 2004 style) and the single array with two tag bits. Every layout
/// must answer "first l-index of an interval" and "next l-index" exactly as the lcp-interval tree of the definition, for every interval.
/// </summary>
[TestFixture]
public class ChildTableTests
{
    private static readonly ChildLayoutKind[] Layouts = Enum.GetValues<ChildLayoutKind>();

    private static ChildTable Build(string text, ChildLayoutKind layout, out int[] lcp)
    {
        int[] sa = NaiveSuffixArray.Build(text);
        lcp = NaiveLcp.Build(text, sa);
        return ChildTable.Build(layout, text.Length, CompactLcp.FromDense(lcp), out _);
    }

    private static void AssertChildrenMatchDefinition(string text, ChildLayoutKind layout, string context)
    {
        ChildTable table = Build(text, layout, out int[] lcp);
        foreach (BruteInterval x in BruteLcpTree.Intervals(lcp))
        {
            int first = table.FirstLIndex(x.Lb, x.Rb, x.LastChild);
            Assert.That(first, Is.EqualTo(x.LIndices[0]), $"{context} {layout} first of [{x.Lb}..{x.Rb}]");
            var chain = new List<int> { first };
            for (int k = table.NextLIndex(first, x.Rb); k >= 0; k = table.NextLIndex(k, x.Rb))
            {
                chain.Add(k);
                Assert.That(chain.Count, Is.LessThanOrEqualTo(x.LIndices.Count + 1), $"{context} {layout} chain of [{x.Lb}..{x.Rb}] does not end");
            }
            Assert.That(chain, Is.EqualTo(x.LIndices), $"{context} {layout} chain of [{x.Lb}..{x.Rb}]");
            Assert.That(lcp[first], Is.EqualTo(x.Ell), $"{context} {layout} ell of [{x.Lb}..{x.Rb}]");
        }
    }

    [Test]
    public void Banana_AllLayouts()
    {
        foreach (ChildLayoutKind layout in Layouts)
        {
            AssertChildrenMatchDefinition("banana", layout, "banana");
            AssertChildrenMatchDefinition("mississippi", layout, "mississippi");
        }
    }

    [Test]
    public void TinyTexts_AllLayouts()
    {
        foreach (ChildLayoutKind layout in Layouts)
            foreach (string s in new[] { "a", "aa", "ab", "ba", "aaa", "aba", "abab" })
                AssertChildrenMatchDefinition(s, layout, $"'{s}'");
    }

    [Test]
    public void AllSmallStrings_AllLayouts()
    {
        foreach (ChildLayoutKind layout in Layouts)
            foreach (var (sigma, maxLen) in new[] { (2, 9), (3, 6), (4, 4) })
                foreach (string s in ExhaustiveStrings.All(sigma, maxLen).Where(s => s.Length > 0))
                    AssertChildrenMatchDefinition(s, layout, $"'{s}'");
    }

    [Test]
    public void AdversarialAndRandomTexts_AllLayouts()
    {
        var rng = new SplitMix64(0x1B02);
        foreach (ChildLayoutKind layout in Layouts)
            for (int i = 0; i < 25; i++)
            {
                int n = 1 + rng.NextInt(300);
                ulong seed = rng.Next();
                foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 3 == 0))
                    AssertChildrenMatchDefinition(text, layout, $"{name} n {n}");
            }
    }

    [Test]
    public void LayoutsDifferOnlyInBytes_ThePackedLayoutsUseFourBytesPerEntry()
    {
        string text = StringFamilies.RandomDna(1000, 9);
        foreach (ChildLayoutKind layout in Layouts)
        {
            ChildTable table = Build(text, layout, out _);
            long expected = layout == ChildLayoutKind.Reference3 ? 12L * 1001 : 4L * 1001;
            Assert.That(table.PayloadBytes, Is.EqualTo(expected), layout.ToString());
        }
    }

    [Test]
    public void TaggedLayout_RefusesTextsThatDoNotFitThirtyBits()
    {
        Assert.That(TaggedChild.MaxLength, Is.EqualTo((1 << 30) - 1));
        Assert.That(() => TaggedChild.Create(TaggedChild.MaxLength + 1), Throws.TypeOf<NotSupportedException>());
    }

    [Test]
    public void TaggedFields_KeepAll30ValueBitsAndTheTagApart()
    {
        // positions near 2^29 and 2^30 cannot be reached by a real text in a unit test, so the field encoding is checked on its own
        foreach (int v in new[] { 1, 2, (1 << 29) - 1, 1 << 29, (1 << 29) + 1, (1 << 30) - 2, TaggedChild.MaxLength })
        {
            Assert.That(TaggedChild.ValueOf(TaggedChild.PackNext(v)), Is.EqualTo(v), $"next {v}");
            Assert.That(TaggedChild.ValueOf(TaggedChild.PackUp(v)), Is.EqualTo(v), $"up {v}");
            Assert.That(TaggedChild.ValueOf(TaggedChild.PackDown(v)), Is.EqualTo(v), $"down {v}");
            Assert.That(TaggedChild.IsNextLink(TaggedChild.PackNext(v)), Is.True, $"next {v}");
            Assert.That(TaggedChild.IsNextLink(TaggedChild.PackUp(v)), Is.False, $"up {v}");
            Assert.That(TaggedChild.IsNextLink(TaggedChild.PackDown(v)), Is.False, $"down {v}");
        }
        Assert.That(TaggedChild.IsNextLink(0), Is.False);
        Assert.That(TaggedChild.PackUp(5), Is.Not.EqualTo(TaggedChild.PackDown(5)));
        Assert.That(TaggedChild.PackNext(5), Is.Not.EqualTo(TaggedChild.PackUp(5)));
    }

    [Test]
    public void RootOfASingleCharacterText_HasItsTwoChildrenAndNothingElse()
    {
        foreach (ChildLayoutKind layout in Layouts)
        {
            ChildTable table = Build("a", layout, out _);   // SA = 1 0, LCP = 0 0: root [0..1] with l-index 1
            Assert.That(table.FirstLIndex(0, 1, true), Is.EqualTo(1), layout.ToString());
            Assert.That(table.NextLIndex(1, 1), Is.EqualTo(-1), layout.ToString());
        }
    }
}
