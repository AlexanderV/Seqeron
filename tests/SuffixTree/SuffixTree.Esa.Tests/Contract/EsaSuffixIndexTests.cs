using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Contract;

/// <summary>
/// The enhanced suffix array behind the index contract (work package 1b.5): <see cref="EsaSuffixIndex"/>. Every operation on every kind of locus is compared with
/// <see cref="IndexContractOracle"/> (derived from the definitions) by <see cref="EsaContractChecker"/>; the tests here add the hand-derived tree of "banana", the capability
/// table of ADR-001 section 6, the rule that core operations never build the navigation pack, argument validation, the three child-table layouts and concurrent first use.
/// </summary>
[TestFixture]
public class EsaSuffixIndexTests
{
    internal static readonly ChildLayoutKind[] Layouts = [ChildLayoutKind.Tagged, ChildLayoutKind.Packed, ChildLayoutKind.Reference3];

    internal static readonly string[] SmallTexts =
    [
        "", "a", "aa", "ab", "ba", "aaa", "aba", "abab", "banana", "mississippi", "abracadabra", "aaaaaaaa", "abababab", "abaababaabaab", "xyzxyz", "abcabcabc",
        "\uD800a\uD800", "\0￿\0", "ab\0￿ab\0￿", "￿￿\0", "mmmmmmm", "abcdefg", "gfedcba", "aabaabaab",
    ];

    internal static (EsaSuffixIndex Index, EsaIndex Esa) Make(string text, ChildLayoutKind layout = ChildLayoutKind.Tagged)
    {
        EsaIndex esa = EsaIndex.Build(text, layout);
        return (new EsaSuffixIndex(esa), esa);
    }

    private static IndexLocus Reach(EsaSuffixIndex index, string pattern)
    {
        Assert.That(index.TryExtendRight(index.Root, pattern.AsSpan(), out IndexLocus locus), Is.True, pattern);
        return locus;
    }

    // ---- construction ------------------------------------------------------------------------------------------------------------

    [Test]
    public void Construction_RejectsNull()
    {
        Assert.That(() => new EsaSuffixIndex(null!), Throws.TypeOf<ArgumentNullException>());
        Assert.That(() => EsaSuffixIndex.Build(null!), Throws.TypeOf<ArgumentNullException>());
    }

    [Test]
    public void Metadata_EngineNameTextAndRoot()
    {
        var (index, esa) = Make("banana");
        Assert.That(index.Engine, Is.EqualTo("esa"));
        Assert.That(index.TextLength, Is.EqualTo(6));
        Assert.That(index.Text.Length, Is.EqualTo(6));
        IIndexText firstView = index.Text;
        Assert.That(index.Text, Is.SameAs(firstView), "one text view per index");
        Assert.That(index.Root, Is.EqualTo(new IndexLocus(0, 6, 0, 0)));
        Assert.That(index.Count(index.Root), Is.EqualTo(6));
        Assert.That(esa.NavigationBuilt, Is.False);
        Assert.That(EsaSuffixIndex.Build("banana", ChildLayoutKind.Packed).TextLength, Is.EqualTo(6));
    }

    [Test]
    public void Text_CopyToAndSymbolAtFollowTheContract()
    {
        var (index, _) = Make("mississippi");
        var buffer = new char[20];
        Assert.That(index.Text.CopyTo(0, buffer), Is.EqualTo(11));
        Assert.That(new string(buffer, 0, 11), Is.EqualTo("mississippi"));
        Assert.That(index.Text.CopyTo(9, buffer.AsSpan(0, 5)), Is.EqualTo(2));
        Assert.That(index.Text.CopyTo(11, buffer), Is.EqualTo(0));
        Assert.That(index.Text.CopyTo(3, Span<char>.Empty), Is.EqualTo(0));
        Assert.That(index.Text.CopyTo(4, buffer.AsSpan(0, 3)), Is.EqualTo(3));
        Assert.That(new string(buffer, 0, 3), Is.EqualTo("iss"));
        Assert.That(() => index.Text.CopyTo(12, buffer), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.Text.CopyTo(-1, buffer), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.Text.SymbolAt(-1), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(index.Text.SymbolAt(0), Is.EqualTo('m'));
        Assert.That(index.Text.SymbolAt(10), Is.EqualTo('i'));
        Assert.That(index.Text.SymbolAt(11), Is.EqualTo(-1));
        Assert.That(index.Text.SymbolAt(long.MaxValue), Is.EqualTo(-1));
    }

    // ---- hand-derived tree of "banana" -------------------------------------------------------------------------------------------

    // text$ = banana$; SA = [6 5 3 1 0 4 2], LCP = [0 0 1 3 0 0 2]; suffixes by rank: $ a ana anana banana na nana.
    // Loci are (Lo, Hi, Length, Aux): the lcp-interval ranks, the matched length, the depth of the node at or below.
    [Test]
    public void Banana_HandDerivedArraysAndTree()
    {
        var (index, esa) = Make("banana");
        long[] sa = [6, 5, 3, 1, 0, 4, 2];
        long[] lcp = [0, 0, 1, 3, 0, 0, 2];
        for (int r = 0; r <= 6; r++)
        {
            Assert.That(index.SuffixAt(r), Is.EqualTo(sa[r]), $"SA[{r}]");
            Assert.That(index.LcpAt(r), Is.EqualTo(lcp[r]), $"LCP[{r}]");
            Assert.That(index.RankOf(sa[r]), Is.EqualTo(r), $"ISA[SA[{r}]]");
        }

        IndexLocus root = index.Root;
        var rootKids = new IndexLocus[4];
        Assert.That(index.GetChildren(root, rootKids), Is.EqualTo(4));
        Assert.That(rootKids, Is.EqualTo(new[] { new IndexLocus(0, 0, 1, 1), new IndexLocus(1, 3, 1, 1), new IndexLocus(4, 4, 7, 7), new IndexLocus(5, 6, 2, 2) }));

        var a = new IndexLocus(1, 3, 1, 1);
        var ana = new IndexLocus(2, 3, 3, 3);
        var na = new IndexLocus(5, 6, 2, 2);
        Assert.That(index.TryChild(root, 'a', out IndexLocus gotA), Is.True);
        Assert.That(gotA, Is.EqualTo(a));
        Assert.That(index.TryChild(root, 'n', out IndexLocus gotNa), Is.True);
        Assert.That(gotNa, Is.EqualTo(na));
        Assert.That(index.TryChild(root, -1, out IndexLocus dollar), Is.True);
        Assert.That(dollar, Is.EqualTo(new IndexLocus(0, 0, 1, 1)), "the leaf of the empty suffix");
        Assert.That(index.Count(dollar), Is.EqualTo(0), "the empty suffix is no text position");
        Assert.That(index.TryChild(root, 'x', out _), Is.False);

        var aKids = new IndexLocus[2];
        Assert.That(index.GetChildren(a, aKids), Is.EqualTo(2));
        Assert.That(aKids, Is.EqualTo(new[] { new IndexLocus(1, 1, 2, 2), ana }));
        var anaKids = new IndexLocus[2];
        Assert.That(index.GetChildren(ana, anaKids), Is.EqualTo(2));
        Assert.That(anaKids, Is.EqualTo(new[] { new IndexLocus(2, 2, 4, 4), new IndexLocus(3, 3, 6, 6) }));
        var naKids = new IndexLocus[2];
        Assert.That(index.GetChildren(na, naKids), Is.EqualTo(2));
        Assert.That(naKids, Is.EqualTo(new[] { new IndexLocus(5, 5, 3, 3), new IndexLocus(6, 6, 5, 5) }));

        Assert.That(Reach(index, "a"), Is.EqualTo(a));
        Assert.That(Reach(index, "ana"), Is.EqualTo(ana));
        Assert.That(Reach(index, "na"), Is.EqualTo(na));
        Assert.That(Reach(index, "an"), Is.EqualTo(new IndexLocus(2, 3, 2, 3)));
        Assert.That(Reach(index, "b"), Is.EqualTo(new IndexLocus(4, 4, 1, 7)));
        Assert.That(Reach(index, "ban"), Is.EqualTo(new IndexLocus(4, 4, 3, 7)));
        Assert.That(Reach(index, "banana"), Is.EqualTo(new IndexLocus(4, 4, 6, 7)));
        Assert.That(Reach(index, "nana"), Is.EqualTo(new IndexLocus(6, 6, 4, 5)));
        Assert.That(index.IsExplicitNode(a), Is.True);
        Assert.That(index.IsExplicitNode(Reach(index, "an")), Is.False);
        Assert.That(index.IsExplicitNode(Reach(index, "banana")), Is.False, "the leaf itself carries the terminator");
        Assert.That(index.IsExplicitNode(index.Root), Is.True);

        Assert.That(index.TryParent(ana, out IndexLocus p), Is.True);
        Assert.That(p, Is.EqualTo(a));
        Assert.That(index.TryParent(a, out p), Is.True);
        Assert.That(p, Is.EqualTo(root));
        Assert.That(index.TryParent(Reach(index, "an"), out p), Is.True);
        Assert.That(p, Is.EqualTo(a));
        Assert.That(index.TryParent(Reach(index, "banana"), out p), Is.True);
        Assert.That(p, Is.EqualTo(root));
        Assert.That(index.TryParent(root, out p), Is.False);
        Assert.That(p, Is.EqualTo(default(IndexLocus)));

        Assert.That(index.TrySuffixLink(a, out IndexLocus link), Is.True);
        Assert.That(link, Is.EqualTo(root));
        Assert.That(index.TrySuffixLink(ana, out link), Is.True);
        Assert.That(link, Is.EqualTo(na));
        Assert.That(index.TrySuffixLink(na, out link), Is.True);
        Assert.That(link, Is.EqualTo(a));
        Assert.That(index.TrySuffixLink(root, out link), Is.False);
        Assert.That(index.TrySuffixLink(new IndexLocus(4, 4, 7, 7), out link), Is.False, "a leaf has no suffix link");
        Assert.That(link, Is.EqualTo(default(IndexLocus)));

        Assert.That(index.Lca(2, 3), Is.EqualTo(ana));
        Assert.That(index.Lca(1, 3), Is.EqualTo(a));
        Assert.That(index.Lca(5, 6), Is.EqualTo(na));
        Assert.That(index.Lca(4, 5), Is.EqualTo(root));
        Assert.That(index.Lca(0, 6), Is.EqualTo(root));
        Assert.That(index.Lca(2, 2), Is.EqualTo(new IndexLocus(2, 2, 4, 4)));
        Assert.That(index.Lca(0, 0), Is.EqualTo(new IndexLocus(0, 0, 1, 1)));
        Assert.That(index.Lcp(2, 3), Is.EqualTo(3));
        Assert.That(index.Lcp(3, 2), Is.EqualTo(3));
        Assert.That(index.Lcp(1, 3), Is.EqualTo(1));
        Assert.That(index.Lcp(5, 6), Is.EqualTo(2));
        Assert.That(index.Lcp(0, 6), Is.EqualTo(0));
        Assert.That(index.Lcp(4, 4), Is.EqualTo(6), "the length of the suffix");
        Assert.That(index.Lcp(0, 0), Is.EqualTo(0));

        Assert.That(index.TryExtendLeft(Reach(index, "an"), 'b', out IndexLocus ban), Is.True);
        Assert.That(ban, Is.EqualTo(new IndexLocus(4, 4, 3, 7)));
        Assert.That(index.TryExtendLeft(Reach(index, "an"), 'a', out _), Is.False);
        Assert.That(index.TryExtendLeft(a, 'n', out IndexLocus nA), Is.True);
        Assert.That(nA, Is.EqualTo(na));
        Assert.That(index.TryExtendLeft(ana, 'n', out IndexLocus nAna), Is.True);
        Assert.That(nAna, Is.EqualTo(new IndexLocus(6, 6, 4, 5)));
        Assert.That(index.TryExtendLeft(root, 'a', out IndexLocus viaRoot), Is.True);
        Assert.That(viaRoot, Is.EqualTo(a));

        var positions = new long[4];
        Assert.That(index.Locate(ana, 0, positions), Is.EqualTo(2));
        Assert.That(positions.Take(2), Is.EqualTo(new long[] { 3, 1 }), "suffix-array order");
        Assert.That(index.Count(a), Is.EqualTo(3));
        Assert.That(index.Count(Reach(index, "an")), Is.EqualTo(2));
        Assert.That(esa.NavigationBuilt, Is.True, "TryParent / TrySuffixLink / RankOf / Lcp / Lca / ExtendLeft have been used");
    }

    [Test]
    public void EdgeSymbol_ReadsTheIncomingEdge()
    {
        var (index, _) = Make("banana");
        var ana = new IndexLocus(2, 3, 3, 3);
        Assert.That(index.EdgeSymbol(ana, 0), Is.EqualTo('n'), "edge a -> ana is 'na'");
        Assert.That(index.EdgeSymbol(ana, 1), Is.EqualTo('a'));
        var leaf = new IndexLocus(4, 4, 7, 7);
        Assert.That(Enumerable.Range(0, 7).Select(k => index.EdgeSymbol(leaf, k)), Is.EqualTo("banana".Select(c => (int)c).Append(-1)));
        Assert.That(index.EdgeSymbol(leaf, 7), Is.EqualTo(-1));
        Assert.That(index.EdgeSymbol(index.Root, 0), Is.EqualTo(-1), "the root has no edge");
        Assert.That(index.EdgeSymbol(new IndexLocus(0, 0, 1, 1), 0), Is.EqualTo(-1));
        Assert.That(index.EdgeSymbol(ana, long.MaxValue), Is.EqualTo(-1));
        Assert.That(() => index.EdgeSymbol(ana, -1), Throws.TypeOf<ArgumentOutOfRangeException>());
    }

    // ---- empty text and one symbol ---------------------------------------------------------------------------------------------

    [Test]
    public void EmptyText_HasOnlyTheRoot()
    {
        foreach (var layout in Layouts)
        {
            var (index, _) = Make("", layout);
            Assert.That(index.TextLength, Is.EqualTo(0));
            Assert.That(index.Root, Is.EqualTo(new IndexLocus(0, 0, 0, 0)));
            Assert.That(index.Root, Is.EqualTo(default(IndexLocus)));
            Assert.That(index.Count(index.Root), Is.EqualTo(0));
            Assert.That(index.IsExplicitNode(index.Root), Is.True);
            Assert.That(index.TryParent(index.Root, out _), Is.False);
            Assert.That(index.TryChild(index.Root, -1, out _), Is.False, "the tree of the empty text is the root alone");
            Assert.That(index.GetChildren(index.Root, new IndexLocus[3]), Is.EqualTo(0));
            Assert.That(index.TrySuffixLink(index.Root, out _), Is.False);
            Assert.That(index.EdgeSymbol(index.Root, 0), Is.EqualTo(-1));
            Assert.That(index.TryExtendRight(index.Root, 'a', out _), Is.False);
            Assert.That(index.TryExtendLeft(index.Root, 'a', out _), Is.False);
            Assert.That(index.Locate(index.Root, 0, new long[2]), Is.EqualTo(0));
            Assert.That(index.SuffixAt(0), Is.EqualTo(0));
            Assert.That(index.RankOf(0), Is.EqualTo(0));
            Assert.That(index.LcpAt(0), Is.EqualTo(0));
            Assert.That(index.Lcp(0, 0), Is.EqualTo(0));
            Assert.That(index.Lca(0, 0), Is.EqualTo(index.Root));
            Assert.That(() => index.SuffixAt(1), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => index.RankOf(1), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => index.LcpAt(1), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => index.Lcp(0, 1), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => index.Count(new IndexLocus(0, 0, 1, 1)), Throws.TypeOf<ArgumentException>(), "there is no terminator leaf in the tree of the empty text");
        }
    }

    [Test]
    public void SingleSymbol_RootLeafAndTerminatorLeaf()
    {
        foreach (var layout in Layouts)
        {
            var (index, _) = Make("a", layout);
            Assert.That(index.Root, Is.EqualTo(new IndexLocus(0, 1, 0, 0)));
            Assert.That(index.Count(index.Root), Is.EqualTo(1));
            var kids = new IndexLocus[2];
            Assert.That(index.GetChildren(index.Root, kids), Is.EqualTo(2));
            Assert.That(kids, Is.EqualTo(new[] { new IndexLocus(0, 0, 1, 1), new IndexLocus(1, 1, 2, 2) }));
            Assert.That(index.TryExtendRight(index.Root, 'a', out IndexLocus a), Is.True);
            Assert.That(a, Is.EqualTo(new IndexLocus(1, 1, 1, 2)));
            Assert.That(index.IsExplicitNode(a), Is.False);
            Assert.That(index.TryParent(a, out IndexLocus parent), Is.True);
            Assert.That(parent, Is.EqualTo(index.Root));
            Assert.That(index.TryExtendRight(a, 'a', out _), Is.False);
            Assert.That(index.TryExtendLeft(a, 'a', out _), Is.False);
            Assert.That(index.TryChild(index.Root, 'a', out IndexLocus leaf), Is.True);
            Assert.That(leaf, Is.EqualTo(new IndexLocus(1, 1, 2, 2)));
            Assert.That(index.TryExtendLeft(new IndexLocus(0, 0, 1, 1), 'a', out IndexLocus viaDollar), Is.True, "a$ is the suffix of length one");
            Assert.That(viaDollar, Is.EqualTo(leaf));
            Assert.That(index.TrySuffixLink(leaf, out _), Is.False);
            Assert.That(index.Lca(0, 1), Is.EqualTo(index.Root));
            Assert.That(index.Lca(1, 1), Is.EqualTo(leaf));
            Assert.That(index.Lcp(1, 1), Is.EqualTo(1));
            Assert.That(index.Lcp(0, 1), Is.EqualTo(0));
        }
    }

    // ---- everything against the oracle ----------------------------------------------------------------------------------------------

    [Test]
    public void SmallTexts_EveryOperationOnEveryLocus_AllThreeLayouts()
    {
        foreach (var layout in Layouts)
            foreach (string text in SmallTexts)
            {
                var (index, _) = Make(text, layout);
                EsaContractChecker.CheckAll(index, new IndexContractOracle(text), $"[{layout}] '{text}'", CheckOptions.Thorough);
            }
    }

    [Test]
    public void AdversarialFamilies_N40_EveryOperationOnEveryLocus()
    {
        foreach (var layout in Layouts)
            foreach (var (name, text) in StringFamilies.Adversarial(40, 7).Where((_, i) => i % 2 == (int)layout % 2))
            {
                var (index, _) = Make(text, layout);
                EsaContractChecker.CheckAll(index, new IndexContractOracle(text), $"[{layout}] {name}", CheckOptions.Thorough);
            }
    }

    [Test]
    public void LongTexts_RepetitiveFamilies_SampledLoci()
    {
        // a^n, Fibonacci, Thue-Morse and periodic texts of a few hundred symbols: LCP values above 255 (escapes of the compact LCP), deep trees, many nested intervals
        var texts = new List<(string Name, string Text)>
        {
            ("a^300", StringFamilies.AllSame(300)),
            ("fibonacci-377", StringFamilies.Fibonacci(377)),
            ("thue-morse-300", StringFamilies.ThueMorse(300)),
            ("periodic-7-280", StringFamilies.Periodic(280, 7, 3)),
            ("ab-repeat-300", StringFamilies.AbRepeat(300)),
        };
        foreach (var (name, text) in texts)
        {
            var (index, _) = Make(text);
            EsaContractChecker.CheckAll(index, new IndexContractOracle(text), name, CheckOptions.Sampled);
        }
    }

    [Test]
    public void LongTexts_RandomFamilies_SampledLoci()
    {
        var texts = new List<(string Name, string Text)>
        {
            ("random-dna-250", StringFamilies.RandomDna(250, 11)),
            ("random-protein-200", StringFamilies.RandomProtein(200, 12)),
            ("random-binary-250", StringFamilies.Random(250, 2, 13)),
            ("skewed-220", StringFamilies.Skewed(220, 14)),
            ("bmp-150", StringFamilies.Bmp(150, 15)),
            ("surrogates-150", StringFamilies.SurrogatePairs(150, 16)),
        };
        foreach (var layout in new[] { ChildLayoutKind.Packed, ChildLayoutKind.Reference3 })
            foreach (var (name, text) in texts)
            {
                var (index, _) = Make(text, layout);
                EsaContractChecker.CheckAll(index, new IndexContractOracle(text), $"[{layout}] {name}", CheckOptions.Sampled);
            }
    }

    [Test]
    public void Layouts_GiveTheSameLociForTheSameQueries()
    {
        foreach (string text in new[] { "mississippi", "abracadabraabracadabra", StringFamilies.Fibonacci(100) })
        {
            var indexes = Layouts.Select(l => Make(text, l).Index).ToArray();
            var oracle = new IndexContractOracle(text);
            foreach (int[] w in oracle.Nodes().Where(w => w.Length == 0 || w[^1] >= 0))
            {
                string pattern = string.Concat(w.Select(s => (char)s));
                IndexLocus expected = oracle.LocusOf(w);
                foreach (var index in indexes)
                {
                    Assert.That(index.TryExtendRight(index.Root, pattern.AsSpan(), out IndexLocus got), Is.True);
                    Assert.That(got, Is.EqualTo(expected), IndexContractOracle.Render(w));
                }
            }
        }
    }

    // ---- capability table ----------------------------------------------------------------------------------------------------------

    private static readonly IndexOperation[] PackDependent =
        [IndexOperation.Parent, IndexOperation.SuffixLink, IndexOperation.InverseSuffixArrayAccess, IndexOperation.Lcp, IndexOperation.Lca];

    /// <summary>The ESA column of the capability table of ADR-001 section 6 (the state once the navigation pack exists).</summary>
    private static IndexSupport Adr001Table(IndexOperation op) => op switch
    {
        IndexOperation.Count or IndexOperation.ExtendRight or IndexOperation.LazyLocate or IndexOperation.TextAccess => IndexSupport.Native,
        IndexOperation.ExtendLeft => IndexSupport.Degraded,
        IndexOperation.Child or IndexOperation.Parent or IndexOperation.StringDepth => IndexSupport.Native,
        IndexOperation.SuffixLink => IndexSupport.Native, // "D until the pack exists, then N"
        IndexOperation.SuffixArrayAccess or IndexOperation.InverseSuffixArrayAccess or IndexOperation.LcpAccess => IndexSupport.Native, // ISA: "N (pack)"
        IndexOperation.Lcp or IndexOperation.Lca => IndexSupport.Native, // "D without the RMQ helper, N with it"
        _ => throw new ArgumentOutOfRangeException(nameof(op)),
    };

    [Test]
    public void CapabilityMatrix_EqualsAdr001Section6_OncePackExists()
    {
        var (index, esa) = Make("banana");
        _ = esa.Navigation;
        var all = Enum.GetValues<IndexOperation>();
        Assert.That(all, Has.Length.EqualTo(14));
        foreach (IndexOperation op in all)
        {
            Assert.That(index.Support(op), Is.EqualTo(Adr001Table(op)), op.ToString());
            if (index.Support(op) == IndexSupport.Degraded) Assert.That(index.SupportNote(op), Is.Not.Null.And.Not.Empty, op.ToString());
        }
        Assert.That(all.Count(op => index.Support(op) == IndexSupport.Unsupported), Is.EqualTo(0), "the ESA does everything");
    }

    [Test]
    public void CapabilityMatrix_BeforeThePack_OnlyThePackOperationsAreDegraded()
    {
        var (index, esa) = Make("banana");
        Assert.That(esa.NavigationBuilt, Is.False);
        foreach (IndexOperation op in Enum.GetValues<IndexOperation>())
        {
            IndexSupport expected = PackDependent.Contains(op) ? IndexSupport.Degraded : Adr001Table(op);
            Assert.That(index.Support(op), Is.EqualTo(expected), op.ToString());
            if (expected == IndexSupport.Degraded)
            {
                string? note = index.SupportNote(op);
                Assert.That(note, Is.Not.Null.And.Not.Empty, op.ToString());
                Assert.That(note, Does.Contain("O("), $"{op}: the note states the cost");
            }
            else Assert.That(index.SupportNote(op), Is.Null, $"{op}: a native operation needs no cost note");
        }
        Assert.That(esa.NavigationBuilt, Is.False, "asking for the capabilities builds nothing");
    }

    [Test]
    public void CapabilityMatrix_SuffixLinkFlipsFromDegradedToNativeWhenThePackIsBuilt()
    {
        var (index, esa) = Make("mississippi");
        Assert.That(index.Support(IndexOperation.SuffixLink), Is.EqualTo(esa.SuffixLinkSupport));
        Assert.That(index.Support(IndexOperation.SuffixLink), Is.EqualTo(IndexSupport.Degraded));
        Assert.That(index.SupportNote(IndexOperation.SuffixLink), Does.Contain("ISA").And.Contain("O(n)"));
        Assert.That(index.TrySuffixLink(Reach(index, "ssi"), out _), Is.True);
        Assert.That(index.Support(IndexOperation.SuffixLink), Is.EqualTo(IndexSupport.Native));
        Assert.That(index.SupportNote(IndexOperation.SuffixLink), Is.Null);
        foreach (IndexOperation op in PackDependent) Assert.That(index.Support(op), Is.EqualTo(IndexSupport.Native), op.ToString());
        Assert.That(index.Support(IndexOperation.ExtendLeft), Is.EqualTo(IndexSupport.Degraded), "O(width) whatever the state");
        Assert.That(index.SupportNote(IndexOperation.ExtendLeft), Does.Contain("width"));
    }

    [Test]
    public void Support_UndefinedOperationThrows()
    {
        var (index, _) = Make("banana");
        Assert.That(() => index.Support((IndexOperation)14), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.Support((IndexOperation)(-1)), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.SupportNote((IndexOperation)14), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.SupportNote((IndexOperation)(-1)), Throws.TypeOf<ArgumentOutOfRangeException>());
    }

    // ---- the navigation pack is lazy ----------------------------------------------------------------------------------------------

    [Test]
    public void CoreOperations_NeverBuildTheNavigationPack()
    {
        foreach (var layout in Layouts)
        {
            var (index, esa) = Make("abracadabraabracadabra", layout);
            IndexLocus abra = Reach(index, "abra");
            Assert.That(index.Count(abra), Is.EqualTo(4));
            Assert.That(index.TryExtendRight(abra, 'c', out IndexLocus abrac), Is.True);
            Assert.That(index.TryExtendRight(index.Root, "cadabra".AsSpan(), out _), Is.True);
            Assert.That(index.Locate(abra, 0, new long[8]), Is.EqualTo(4));
            Assert.That(index.Locate(abra, 1, new long[2]), Is.EqualTo(2));
            Assert.That(index.Text.SymbolAt(4), Is.EqualTo('c'));
            Assert.That(index.Text.CopyTo(0, new char[5]), Is.EqualTo(5));
            Assert.That(index.IsExplicitNode(abra), Is.True);
            Assert.That(index.IsExplicitNode(abrac), Is.False);
            Assert.That(index.TryChild(abra, 'c', out _), Is.True);
            Assert.That(index.TryChild(abra, -1, out _), Is.True, "the last occurrence of abra ends the text: a terminator child");
            Assert.That(index.TryChild(abra, 'z', out _), Is.False);
            var kids = new IndexLocus[8];
            Assert.That(index.GetChildren(index.Root, kids), Is.GreaterThan(3));
            Assert.That(index.EdgeSymbol(abra, 1), Is.EqualTo('r'), "the edge from the node a is bra");
            Assert.That(index.SuffixAt(3), Is.GreaterThanOrEqualTo(0));
            Assert.That(index.LcpAt(5), Is.GreaterThanOrEqualTo(0));
            Assert.That(index.TryExtendLeft(index.Root, 'a', out _), Is.True, "the empty string extended to the left is the symbol: no ISA needed");
            Assert.That(esa.NavigationBuilt, Is.False, layout.ToString());
            Assert.That(esa.NavigationBuildCount, Is.EqualTo(0), layout.ToString());
            Assert.That(index.Support(IndexOperation.SuffixLink), Is.EqualTo(IndexSupport.Degraded));
        }
    }

    [Test]
    public void EachPackOperation_BuildsThePackExactlyOnce()
    {
        var calls = new Dictionary<string, Action<EsaSuffixIndex>>
        {
            ["TryParent"] = i => Assert.That(i.TryParent(Reach(i, "ana"), out _), Is.True),
            ["TrySuffixLink"] = i => Assert.That(i.TrySuffixLink(Reach(i, "ana"), out _), Is.True),
            ["TryExtendLeft"] = i => Assert.That(i.TryExtendLeft(Reach(i, "ana"), 'b', out _), Is.True),
            ["RankOf"] = i => Assert.That(i.RankOf(2), Is.GreaterThan(0)),
            ["Lcp"] = i => Assert.That(i.Lcp(2, 3), Is.EqualTo(3)),
            ["Lca"] = i => Assert.That(i.Lca(2, 3), Is.Not.EqualTo(i.Root)),
        };
        foreach (var (name, call) in calls)
        {
            var (index, esa) = Make("banana");
            Assert.That(esa.NavigationBuilt, Is.False, name);
            call(index);
            Assert.That(esa.NavigationBuilt, Is.True, name);
            Assert.That(esa.NavigationBuildCount, Is.EqualTo(1), name);
            call(index);
            Assert.That(esa.NavigationBuildCount, Is.EqualTo(1), $"{name}: the pack is built once");
        }
    }

    [Test]
    public void InvalidArguments_NeverBuildThePack()
    {
        var (index, esa) = Make("banana");
        IndexLocus foreign = new(9, 9, 9, 9);
        IndexLocus ana = Reach(index, "ana");
        Assert.That(() => index.TryParent(foreign, out _), Throws.TypeOf<ArgumentException>());
        Assert.That(() => index.TrySuffixLink(foreign, out _), Throws.TypeOf<ArgumentException>());
        Assert.That(() => index.TrySuffixLink(Reach(index, "an"), out _), Throws.TypeOf<ArgumentException>(), "not a node");
        Assert.That(() => index.TryExtendLeft(foreign, 'a', out _), Throws.TypeOf<ArgumentException>());
        Assert.That(index.TryExtendLeft(ana, -1, out _), Is.False);
        Assert.That(index.TryExtendLeft(ana, 65536, out _), Is.False);
        Assert.That(() => index.RankOf(-1), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.RankOf(7), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.RankOf(long.MaxValue), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.Lcp(-1, 2), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.Lcp(2, 7), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.Lca(7, 2), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.Lca(2, -1), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(esa.NavigationBuilt, Is.False);
    }

    // ---- argument validation ----------------------------------------------------------------------------------------------------------

    [Test]
    public void MalformedLoci_AreRejectedByEveryOperation()
    {
        var (index, _) = Make("banana");
        var bad = new[]
        {
            default(IndexLocus),                     // not the root of a non-empty text
            new IndexLocus(0, 7, 0, 0),              // Hi beyond n
            new IndexLocus(-1, 3, 1, 1),             // negative Lo
            new IndexLocus(3, 2, 1, 1),              // Lo > Hi
            new IndexLocus(1, 3, -1, 1),             // negative length
            new IndexLocus(1, 3, 2, 1),              // Aux below the length
            new IndexLocus(1, 1, 1, 8),              // deeper than any suffix with its terminator
            new IndexLocus(1, 6, 0, 0),              // length 0 but not the whole range
            new IndexLocus(0, 3, 2, 2),              // rank 0 inside a non-root locus
            new IndexLocus(0, 0, 2, 2),              // the empty-suffix leaf has depth 1
            new IndexLocus(long.MaxValue, long.MaxValue, 1, 1),
            new IndexLocus(0, long.MaxValue, 0, 0),
        };
        foreach (IndexLocus locus in bad)
        {
            string at = locus.ToString();
            Assert.That(() => index.Count(locus), Throws.TypeOf<ArgumentException>(), at);
            Assert.That(() => index.Locate(locus, 0, new long[2]), Throws.TypeOf<ArgumentException>(), at);
            Assert.That(() => index.TryExtendRight(locus, 'a', out _), Throws.TypeOf<ArgumentException>(), at);
            Assert.That(() => index.TryExtendRight(locus, "a".AsSpan(), out _), Throws.TypeOf<ArgumentException>(), at);
            Assert.That(() => index.TryExtendRight(locus, ReadOnlySpan<char>.Empty, out _), Throws.TypeOf<ArgumentException>(), at);
            Assert.That(() => index.TryExtendLeft(locus, 'a', out _), Throws.TypeOf<ArgumentException>(), at);
            Assert.That(() => index.IsExplicitNode(locus), Throws.TypeOf<ArgumentException>(), at);
            Assert.That(() => index.TryParent(locus, out _), Throws.TypeOf<ArgumentException>(), at);
            Assert.That(() => index.TryChild(locus, 'a', out _), Throws.TypeOf<ArgumentException>(), at);
            Assert.That(() => index.GetChildren(locus, new IndexLocus[4]), Throws.TypeOf<ArgumentException>(), at);
            Assert.That(() => index.TrySuffixLink(locus, out _), Throws.TypeOf<ArgumentException>(), at);
            Assert.That(() => index.EdgeSymbol(locus, 0), Throws.TypeOf<ArgumentException>(), at);
        }
    }

    private static bool ShapeAccepted(int n, long lo, long hi, long length, long aux)
    {
        // the O(1) shape rules documented on EsaSuffixIndex: inside the index, aux >= length, the root, the leaf of the empty suffix or a point of a non-empty suffix
        if (lo < 0 || lo > hi || hi > n || length < 0 || aux < length || aux > n + 1L) return false;
        return length == 0 ? lo == 0 && hi == n && aux == 0 : n > 0 && (lo > 0 || (hi == 0 && aux == 1));
    }

    private static bool NodeAccepted(IndexContractOracle o, EsaIndex esa, int lo, int hi, long length, long aux)
    {
        // the checks of the node operations, from the oracle's LCP array: both borders of the interval below the depth, the matched length below it and above the parent's depth, a leaf has the depth
        // of its suffix with the terminator, an inner interval has its first l-index in the child table (inside the interval) with the LCP value of the depth
        int n = o.N;
        if (aux == 0) return true;
        int left = o.Lcp[lo];
        int right = hi < n ? o.Lcp[hi + 1] : -1;
        if (Math.Max(left, right) >= length) return false;
        if (lo == hi) return aux == n - o.Sa[lo] + 1;
        int first = esa.Child.FirstLIndex(lo, hi, left > right);
        return first > lo && first <= hi && o.Lcp[first] == aux;
    }

    [Test]
    public void ForgedLoci_CountAndTheNodeOperationsAcceptExactlyTheDocumentedTuples()
    {
        // every tuple in a box around the valid loci: Count / Locate accept exactly the tuples of the shape rules, the node operations exactly those that also pass the node checks. Every real locus passes them,
        // and a single-rank tuple (a point of a leaf edge) passes only when it is real: the depth of a leaf is known. For wider intervals the O(1) checks cannot exclude every interval that merely looks like a
        // node (the interior is not scanned), so look-alikes pass: they are counted and reported; a locus that this index hands out is never one.
        foreach (var (text, layout) in new[] { ("banana", ChildLayoutKind.Tagged), ("mississippi", ChildLayoutKind.Packed), ("abababa", ChildLayoutKind.Reference3), ("aaaa", ChildLayoutKind.Tagged), ("a", ChildLayoutKind.Packed), ("", ChildLayoutKind.Reference3) })
        {
            var (index, esa) = Make(text, layout);
            var oracle = new IndexContractOracle(text);
            int n = text.Length;
            var valid = new HashSet<IndexLocus>();
            foreach (int[] node in oracle.Nodes())
            {
                valid.Add(oracle.LocusOf(node));
                int[] parent = oracle.Parent(node) ?? [];
                for (int len = parent.Length + 1; len <= node.Length; len++) valid.Add(oracle.LocusOf(node[..len]));
            }
            int accepted = 0;
            var positions = new long[n + 2];
            for (int lo = -1; lo <= n + 1; lo++)
                for (int hi = -1; hi <= n + 1; hi++)
                    for (long length = -1; length <= n + 2; length++)
                        for (long aux = -1; aux <= n + 2; aux++)
                        {
                            var locus = new IndexLocus(lo, hi, length, aux);
                            string at = $"[{layout}] '{text}' {locus}";
                            bool shape = ShapeAccepted(n, lo, hi, length, aux);
                            bool nodeOk = shape && NodeAccepted(oracle, esa, lo, hi, length, aux);
                            if (valid.Contains(locus)) Assert.That(nodeOk, Is.True, at + " is a locus");
                            if (shape)
                            {
                                Assert.That(index.Count(locus), Is.EqualTo(hi - Math.Max(lo, 1) + 1), at);
                                Assert.That(index.Locate(locus, 0, positions), Is.EqualTo(hi - Math.Max(lo, 1) + 1), at);
                            }
                            else
                            {
                                Assert.That(() => index.Count(locus), Throws.TypeOf<ArgumentException>(), at);
                                Assert.That(() => index.Locate(locus, 0, positions), Throws.TypeOf<ArgumentException>(), at);
                            }
                            if (nodeOk)
                            {
                                accepted++;
                                Assert.That(index.IsExplicitNode(locus), Is.EqualTo(length == aux), at);
                                Assert.That(index.TryExtendRight(locus, ReadOnlySpan<char>.Empty, out IndexLocus same), Is.True, at);
                                Assert.That(same, Is.EqualTo(locus), at);
                                if (lo == hi) Assert.That(valid.Contains(locus), Is.True, at + " is no locus of the text");
                            }
                            else
                            {
                                Assert.That(() => index.IsExplicitNode(locus), Throws.TypeOf<ArgumentException>(), at);
                                Assert.That(() => index.TryExtendRight(locus, ReadOnlySpan<char>.Empty, out _), Throws.TypeOf<ArgumentException>(), at);
                            }
                        }
            Assert.That(accepted, Is.GreaterThanOrEqualTo(valid.Count), text);
            TestContext.Out.WriteLine($"[{layout}] '{text}': {valid.Count} loci, {accepted - valid.Count} look-alikes pass the node checks");
        }
    }

    [Test]
    public void InvalidSymbolWithAMalformedLocus_StillThrows()
    {
        var (index, _) = Make("banana");
        var bad = new IndexLocus(9, 9, 9, 9);
        Assert.That(() => index.TryExtendRight(bad, -1, out _), Throws.TypeOf<ArgumentException>());
        Assert.That(() => index.TryExtendRight(bad, 70000, out _), Throws.TypeOf<ArgumentException>());
        Assert.That(() => index.TryExtendLeft(bad, -1, out _), Throws.TypeOf<ArgumentException>());
        Assert.That(() => index.TryExtendLeft(bad, 70000, out _), Throws.TypeOf<ArgumentException>());
    }

    // ---- corrupt child tables end in InvalidOperationException, never in a loop or a read outside the interval ---------------------------------

    private static T Field<T>(object owner, string name)
        => (T)owner.GetType().GetField(name, System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)!.GetValue(owner)!;

    private static void AssertWalksRefuse(EsaSuffixIndex index, string what)
    {
        Assert.Throws<InvalidOperationException>(() => index.GetChildren(index.Root, new IndexLocus[8]), $"{what}: GetChildren");
        Assert.Throws<InvalidOperationException>(() => index.GetChildren(index.Root, Span<IndexLocus>.Empty), $"{what}: GetChildren (count only)");
        Assert.Throws<InvalidOperationException>(() => index.TryChild(index.Root, 'z', out _), $"{what}: TryChild of a key behind every child");
        Assert.Throws<InvalidOperationException>(() => index.TryChild(index.Root, 'b', out _), $"{what}: TryChild of a key in the middle");
    }

    [Test]
    public void CorruptChildTable_Reference3_IsRefused()
    {
        const string text = "abracadabra";
        EsaIndex esa = EsaIndex.Build(text, ChildLayoutKind.Reference3);
        var index = new EsaSuffixIndex(esa);
        object table = esa.Child.Reference3;
        int[] next = Field<int[]>(table, "_next");
        int[] up = Field<int[]>(table, "_up");
        int[] down = Field<int[]>(table, "_down");
        int first = esa.Child.FirstLIndex(0, esa.Length, true);
        int second = next[first];

        next[first] = first;                       // the l-index after the first points at itself
        AssertWalksRefuse(index, "self loop");
        next[first] = esa.Length + 3;              // ... or beyond the interval
        AssertWalksRefuse(index, "boundary beyond the interval");
        next[first] = second;

        int saved = down[0];
        down[0] = 0;                               // the first l-index of the root is its own left border
        AssertWalksRefuse(index, "first boundary at the left border");
        down[0] = saved;

        var kids = new IndexLocus[6];
        Assert.That(index.GetChildren(index.Root, kids), Is.EqualTo(6), "the repaired table works again");
        IndexLocus a = kids[1]; // the node 'a': a non-last child, its first l-index sits at up[rb]
        int savedUp = up[a.Hi];
        up[a.Hi] = 0;                              // the first l-index of an inner child at its left border
        Assert.Throws<InvalidOperationException>(() => index.GetChildren(index.Root, new IndexLocus[8]), "child with a bad first l-index: GetChildren");
        Assert.Throws<InvalidOperationException>(() => index.TryChild(index.Root, 'a', out _), "child with a bad first l-index: TryChild");
        up[a.Hi] = (int)a.Hi + 5;                       // ... or behind its interval
        Assert.Throws<InvalidOperationException>(() => index.GetChildren(index.Root, new IndexLocus[8]), "child with a first l-index behind its interval");
        up[a.Hi] = savedUp;
        Assert.That(index.GetChildren(index.Root, kids), Is.EqualTo(6));
    }

    [Test]
    public void CorruptChildTable_Tagged_IsRefused()
    {
        const string text = "abracadabra";
        EsaIndex esa = EsaIndex.Build(text, ChildLayoutKind.Tagged);
        var index = new EsaSuffixIndex(esa);
        object table = esa.Child.Tagged;
        uint[] cld = Field<uint[]>(table, "_cld");
        int first = esa.Child.FirstLIndex(0, esa.Length, true);
        uint savedNext = cld[first];
        uint savedRoot = cld[0];

        cld[first] = TaggedChild.PackNext(first);
        AssertWalksRefuse(index, "self loop");
        cld[first] = savedNext;

        cld[0] = TaggedChild.PackDown(0);
        AssertWalksRefuse(index, "first boundary at the left border");
        cld[0] = savedRoot;

        var kids = new IndexLocus[6];
        Assert.That(index.GetChildren(index.Root, kids), Is.EqualTo(6), "the repaired table works again");
        IndexLocus a = kids[1];
        uint savedUp = cld[a.Hi];
        cld[a.Hi] = TaggedChild.PackUp(0);
        Assert.Throws<InvalidOperationException>(() => index.GetChildren(index.Root, new IndexLocus[8]), "child with a bad first l-index: GetChildren");
        Assert.Throws<InvalidOperationException>(() => index.TryChild(index.Root, 'a', out _), "child with a bad first l-index: TryChild");
        cld[a.Hi] = savedUp;
        Assert.That(index.GetChildren(index.Root, kids), Is.EqualTo(6));
    }

    [Test]
    public void NodeOperations_RejectPointsInsideAnEdge()
    {
        var (index, _) = Make("banana");
        IndexLocus mid = Reach(index, "an");
        Assert.That(index.IsExplicitNode(mid), Is.False);
        Assert.That(() => index.TryChild(mid, 'a', out _), Throws.TypeOf<ArgumentException>());
        Assert.That(() => index.GetChildren(mid, new IndexLocus[4]), Throws.TypeOf<ArgumentException>());
        Assert.That(() => index.TrySuffixLink(mid, out _), Throws.TypeOf<ArgumentException>());
        Assert.That(() => index.EdgeSymbol(mid, 0), Throws.TypeOf<ArgumentException>());
        Assert.That(index.TryParent(mid, out _), Is.True, "the parent of a point is defined");
    }

    [Test]
    public void ArrayAccessors_RejectOutOfRangeArguments()
    {
        var (index, _) = Make("banana");
        foreach (long bad in new long[] { -1, 7, long.MinValue, long.MaxValue, 1L << 32, (1L << 32) + 3 })
        {
            Assert.That(() => index.SuffixAt(bad), Throws.TypeOf<ArgumentOutOfRangeException>(), bad.ToString());
            Assert.That(() => index.LcpAt(bad), Throws.TypeOf<ArgumentOutOfRangeException>(), bad.ToString());
            Assert.That(() => index.RankOf(bad), Throws.TypeOf<ArgumentOutOfRangeException>(), bad.ToString());
            Assert.That(() => index.Lcp(bad, 1), Throws.TypeOf<ArgumentOutOfRangeException>(), bad.ToString());
            Assert.That(() => index.Lcp(1, bad), Throws.TypeOf<ArgumentOutOfRangeException>(), bad.ToString());
            Assert.That(() => index.Lca(bad, 1), Throws.TypeOf<ArgumentOutOfRangeException>(), bad.ToString());
            Assert.That(() => index.Lca(1, bad), Throws.TypeOf<ArgumentOutOfRangeException>(), bad.ToString());
        }
        Assert.That(index.SuffixAt(6), Is.EqualTo(2));
        Assert.That(index.RankOf(6), Is.EqualTo(0), "ISA[n] is the empty suffix");
        Assert.That(index.LcpAt(6), Is.EqualTo(2));
    }

    [Test]
    public void Locate_NegativeSkipAndExhaustedPages()
    {
        var (index, _) = Make("banana");
        IndexLocus a = Reach(index, "a");
        Assert.That(() => index.Locate(a, -1, new long[3]), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.Locate(a, -1, Span<long>.Empty), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.Locate(a, long.MinValue, new long[3]), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.Locate(new IndexLocus(0, 0, 1, 1), -1, new long[3]), Throws.TypeOf<ArgumentOutOfRangeException>(), "also for a locus without positions");
        Assert.That(index.Locate(a, 3, new long[3]), Is.EqualTo(0));
        Assert.That(index.Locate(a, long.MaxValue, new long[3]), Is.EqualTo(0));
        Assert.That(index.Locate(a, 2, new long[3]), Is.EqualTo(1));
        Assert.That(index.Locate(a, 0, Span<long>.Empty), Is.EqualTo(0));
        Assert.That(index.Locate(new IndexLocus(0, 0, 1, 1), 0, new long[3]), Is.EqualTo(0), "the leaf of the empty suffix has no text position");
        var all = new long[10];
        Assert.That(index.Locate(index.Root, 0, all), Is.EqualTo(6), "the root locates every text position and not n");
        Assert.That(all.Take(6), Is.EqualTo(new long[] { 5, 3, 1, 0, 4, 2 }));
        Assert.That(index.Locate(index.Root, 4, all), Is.EqualTo(2));
    }

    [Test]
    public void ExtendRight_LongPatternsAndEmptyPatternOnPoints()
    {
        var (index, _) = Make(StringFamilies.Fibonacci(2000));
        string text = StringFamilies.Fibonacci(2000);
        Assert.That(index.TryExtendRight(index.Root, text.AsSpan(), out IndexLocus whole), Is.True);
        Assert.That(whole.Length, Is.EqualTo(2000));
        Assert.That(index.Count(whole), Is.EqualTo(1));
        Assert.That(index.TryExtendRight(index.Root, (text + "a").AsSpan(), out _), Is.False);
        Assert.That(index.TryExtendRight(whole, ReadOnlySpan<char>.Empty, out IndexLocus same), Is.True);
        Assert.That(same, Is.EqualTo(whole));
        Assert.That(index.TryExtendRight(index.Root, text.AsSpan(1500), out IndexLocus tail), Is.True);
        Assert.That(index.Count(tail), Is.GreaterThanOrEqualTo(1));
    }

    // ---- concurrency ---------------------------------------------------------------------------------------------------------------------

    [Test]
    public void ConcurrentFirstUse_BuildsThePackOnceAndEveryReaderGetsTheSameAnswers()
    {
        string text = StringFamilies.RandomDna(400, 5);
        var oracle = new IndexContractOracle(text);
        var (index, esa) = Make(text);
        var probes = oracle.Nodes().Where(w => w.Length > 0 && w[^1] >= 0).Take(120).ToArray();
        IndexLocus[] loci = probes.Select(oracle.LocusOf).ToArray();
        // expected answers from the oracle, never from the index under test
        var expectedLinks = probes.Select(w => oracle.LocusOf(w[1..])).ToArray();
        var expectedParents = probes.Select(w => oracle.LocusOf(oracle.Parent(w)!)).ToArray();
        var errors = new System.Collections.Concurrent.ConcurrentBag<string>();
        using var start = new ManualResetEventSlim(false);
        var threads = Enumerable.Range(0, 8).Select(t => new Thread(() =>
        {
            start.Wait();
            try
            {
                for (int round = 0; round < 3; round++)
                    for (int k = 0; k < loci.Length; k++)
                    {
                        int at = (k + t * 13) % loci.Length;
                        if (!index.TrySuffixLink(loci[at], out IndexLocus link) || link != expectedLinks[at]) errors.Add($"link {at}");
                        if (!index.TryParent(loci[at], out IndexLocus parent) || parent != expectedParents[at]) errors.Add($"parent {at}");
                        if (index.RankOf(at) != oracle.Isa[at]) errors.Add($"rank {at}");
                        if (index.Lcp(at, 200 + at) != oracle.LcpOfRanks(at, 200 + at)) errors.Add($"lcp {at}");
                    }
            }
            catch (Exception e)
            {
                errors.Add(e.ToString());
            }
        })).ToArray();
        foreach (var thread in threads) thread.Start();
        start.Set();
        foreach (var thread in threads) thread.Join();
        Assert.That(errors, Is.Empty);
        Assert.That(esa.NavigationBuildCount, Is.EqualTo(1));
        Assert.That(index.Support(IndexOperation.SuffixLink), Is.EqualTo(IndexSupport.Native));
    }

    [Test]
    public void ConcurrentReadersOfTheCoreAgree()
    {
        string text = StringFamilies.RandomDna(500, 9);
        var oracle = new IndexContractOracle(text);
        var (index, esa) = Make(text);
        var patterns = Enumerable.Range(0, 100).Select(i => text.Substring(i * 3, 2 + (i % 7))).ToArray();
        long[] expected = patterns.Select(p => (long)BruteForce.Occurrences(text, p).Count).ToArray();
        var errors = new System.Collections.Concurrent.ConcurrentBag<string>();
        Parallel.For(0, 6, new ParallelOptions { MaxDegreeOfParallelism = 6 }, _ =>
        {
            for (int k = 0; k < patterns.Length; k++)
            {
                if (!index.TryExtendRight(index.Root, patterns[k].AsSpan(), out IndexLocus locus) || index.Count(locus) != expected[k]) errors.Add(patterns[k]);
                if (!index.IsExplicitNode(locus) && (!index.TryParent(locus, out IndexLocus parent) || index.Count(parent) < expected[k])) errors.Add("parent " + patterns[k]);
            }
        });
        Assert.That(errors, Is.Empty);
        Assert.That(esa.NavigationBuildCount, Is.LessThanOrEqualTo(1));
    }
}
