using System.Runtime.InteropServices;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Contract;

/// <summary>The engines of the aliasing suites.</summary>
public static class AliasEngines
{
    /// <summary>Every engine that answers the navigation members with enhanced-suffix-array loci: the in-memory index with its three child-table layouts, the facade with two of them, and the container reader at both access modes (the Ukkonen adapter of the oracle has no navigation; the property test covers its extension members).</summary>
    public static IEnumerable<IndexEngineFactory> Navigation { get; } = IndexEngineFactories.All.Where(f => f.Name != "ukkonen-adapter").ToArray();
}

/// <summary>
/// The result of a locus-returning member may be the very variable that holds the locus (work package 1b.12; finding R01 of the review of 1b.10, where <c>TryExtendRight</c> returned false for a string that occurs, the span overload
/// returned a locus of the wrong depth and <c>TryExtendLeft</c> threw for every non-root locus, because the implementations cleared <c>result</c> before they read <c>locus</c>). One test per member and overload: for every locus of three small
/// texts (the root, the points inside edges, the explicit nodes, the leaves with and without the terminator, the leaf of the empty suffix) and every argument, the call with ONE variable for the locus and the result, the call with two
/// variables and the answer of the definitional oracle are the same, and so is the failing answer (false, the default locus) and the exception of a locus that is not an explicit node. <c>GetChildren</c> writes its children into the
/// array that holds the node, <c>Locate</c> writes its positions over the bytes of the locus. Parameterised over the engines that have the navigation members.
/// </summary>
[TestFixtureSource(typeof(AliasEngines), nameof(AliasEngines.Navigation))]
public class IndexAliasingTests(IndexEngineFactory factory)
{
    private static readonly string[] Texts = ["abracadabra", "mississippi", "aaaaaaa"];

    private static readonly Dictionary<string, IndexContractOracle> Oracles = [];

    private static IndexContractOracle OracleOf(string text)
    {
        lock (Oracles)
        {
            if (Oracles.TryGetValue(text, out IndexContractOracle? known)) return known;
            var oracle = new IndexContractOracle(text);
            Oracles.Add(text, oracle);
            return oracle;
        }
    }

    private ISuffixTreeNavigation Nav(string text) => (ISuffixTreeNavigation)factory.Create(text);

    private static void Same(Outcome actual, Outcome expected, Func<string> what)
    {
        if (!actual.Equals(expected)) Assert.Fail($"{what()}: expected {expected}, got {actual}");
    }

    private static string At(string text, int[] w, string call) => $"'{text}' {IndexContractOracle.Render(w)} {call}";

    /// <summary>For every text, every locus and every argument: the separate-variable call, the aliased call and the oracle agree.</summary>
    private void ForEveryLocus(Action<string, IndexContractOracle, ISuffixTreeNavigation, int[], IndexLocus> check)
    {
        foreach (string text in Texts)
        {
            IndexContractOracle o = OracleOf(text);
            ISuffixTreeNavigation nav = Nav(text);
            foreach (int[] w in AliasCatalog.Loci(o)) check(text, o, nav, w, o.LocusOf(w));
        }
    }

    [Test]
    public void TheResultMayAliasTheLocus_AdvancingOneVariableAlongAString()
    {
        ISuffixTreeNavigation nav = Nav("abracadabra");
        IndexLocus apart = nav.Root, alias = nav.Root;
        foreach (char c in "abr")
        {
            Assert.That(nav.TryExtendRight(apart, c, out IndexLocus next), Is.True);
            apart = next;
            Assert.That(nav.TryExtendRight(in alias, c, out alias), Is.True, $"aliased step '{c}'");
            Assert.That(alias, Is.EqualTo(apart), $"aliased step '{c}'");
        }
        Assert.That(alias.Length, Is.EqualTo(3));
        Assert.That(nav.Count(alias), Is.EqualTo(2));
        Assert.That(nav.TryExtendLeft(apart, 'd', out IndexLocus leftApart), Is.True, "dabr occurs");
        IndexLocus left = apart;
        Assert.That(nav.TryExtendLeft(in left, 'd', out left), Is.True);
        Assert.That(left, Is.EqualTo(leftApart));
        Assert.That(left.Length, Is.EqualTo(4));
        Assert.That(nav.TryExtendLeft(apart, 'a', out leftApart), Is.False, "aabr does not occur");
        left = apart;
        Assert.That(nav.TryExtendLeft(in left, 'a', out left), Is.False);
        Assert.That(left, Is.EqualTo(default(IndexLocus)).And.EqualTo(leftApart));
    }

    [Test]
    public void TheNaturalLoop_OneVariableAdvancedAlongEverySubstring_ReachesTheLocusOfTheOracle()
    {
        foreach (string text in Texts)
        {
            IndexContractOracle o = OracleOf(text);
            ISuffixTreeNavigation nav = Nav(text);
            foreach (string pattern in AliasCatalog.Patterns(text, 8).Where(p => p.Length > 0 && text.Contains(p, StringComparison.Ordinal)))
            {
                IndexLocus bySymbol = nav.Root, byPattern = nav.Root;
                int length = 0;
                foreach (char c in pattern)
                {
                    length++;
                    int[] prefix = IndexContractOracle.Symbols(pattern[..length]);
                    Assert.That(nav.TryExtendRight(in bySymbol, c, out bySymbol), Is.True, $"'{text}' '{pattern}' symbol {length}");
                    Assert.That(bySymbol, Is.EqualTo(o.LocusOf(prefix)), $"'{text}' '{pattern}' symbol {length}");
                    Assert.That(nav.TryExtendRight(in byPattern, pattern.AsSpan(length - 1, 1), out byPattern), Is.True, $"'{text}' '{pattern}' pattern {length}");
                    Assert.That(byPattern, Is.EqualTo(bySymbol), $"'{text}' '{pattern}' pattern {length}");
                }
                IndexLocus whole = nav.Root;
                Assert.That(nav.TryExtendRight(in whole, pattern.AsSpan(), out whole), Is.True, $"'{text}' '{pattern}' whole");
                Assert.That(whole, Is.EqualTo(bySymbol), $"'{text}' '{pattern}' whole");
            }
        }
    }

    [Test]
    public void TryExtendRight_Symbol_TheResultMayAliasTheLocus()
    {
        ForEveryLocus((text, o, nav, w, locus) =>
        {
            foreach (int c in AliasCatalog.Symbols(o))
            {
                string at = At(text, w, $"TryExtendRight({c})");
                Outcome expected = AliasCatalog.Expect(o, c is >= 0 and <= 0xFFFF ? o.ExtendRight(w, c) : null);
                Same(AliasProbe.ExtendRightApart(nav, locus, c), expected, () => at + " (two variables)");
                Same(AliasProbe.ExtendRightAlias(nav, locus, c), expected, () => at + " (one variable)");
            }
        });
    }

    [Test]
    public void TryExtendRight_Span_TheResultMayAliasTheLocus()
    {
        ForEveryLocus((text, o, nav, w, locus) =>
        {
            foreach (string pattern in AliasCatalog.Patterns(text, 4))
            {
                string at = At(text, w, $"TryExtendRight(\"{pattern}\")");
                int[] extended = [.. w, .. IndexContractOracle.Symbols(pattern)];
                Outcome expected = AliasCatalog.Expect(o, o.Exists(extended) ? extended : null);
                Same(AliasProbe.ExtendRightApart(nav, locus, pattern), expected, () => at + " (two variables)");
                Same(AliasProbe.ExtendRightAlias(nav, locus, pattern), expected, () => at + " (one variable)");
            }
        });
    }

    [Test]
    public void TryExtendLeft_TheResultMayAliasTheLocus()
    {
        ForEveryLocus((text, o, nav, w, locus) =>
        {
            foreach (int c in AliasCatalog.Symbols(o))
            {
                string at = At(text, w, $"TryExtendLeft({c})");
                Outcome expected = AliasCatalog.Expect(o, c is >= 0 and <= 0xFFFF ? o.ExtendLeft(c, w) : null);
                Same(AliasProbe.ExtendLeftApart(nav, locus, c), expected, () => at + " (two variables)");
                Same(AliasProbe.ExtendLeftAlias(nav, locus, c), expected, () => at + " (one variable)");
            }
        });
    }

    [Test]
    public void TryParent_TheResultMayAliasTheLocus()
    {
        ForEveryLocus((text, o, nav, w, locus) =>
        {
            string at = At(text, w, "TryParent");
            Outcome expected = AliasCatalog.Expect(o, o.Parent(w));
            Same(AliasProbe.ParentApart(nav, locus), expected, () => at + " (two variables)");
            Same(AliasProbe.ParentAlias(nav, locus), expected, () => at + " (one variable)");
        });
    }

    [Test]
    public void TryChild_TheResultMayAliasTheNode()
    {
        ForEveryLocus((text, o, nav, w, locus) =>
        {
            foreach (int key in AliasCatalog.Symbols(o))
            {
                string at = At(text, w, $"TryChild({key})");
                if (!o.IsExplicit(w))
                {
                    Assert.That(() => AliasProbe.ChildApart(nav, locus, key), Throws.TypeOf<ArgumentException>(), at + " (two variables) is no node");
                    Assert.That(() => AliasProbe.ChildAlias(nav, locus, key), Throws.TypeOf<ArgumentException>(), at + " (one variable) is no node");
                    continue;
                }
                int[]? child = o.Children(w).Where(k => k.Key == key).Select(k => k.Child).FirstOrDefault();
                Outcome expected = AliasCatalog.Expect(o, child);
                Same(AliasProbe.ChildApart(nav, locus, key), expected, () => at + " (two variables)");
                Same(AliasProbe.ChildAlias(nav, locus, key), expected, () => at + " (one variable)");
            }
        });
    }

    [Test]
    public void TrySuffixLink_TheResultMayAliasTheNode()
    {
        ForEveryLocus((text, o, nav, w, locus) =>
        {
            string at = At(text, w, "TrySuffixLink");
            if (!o.IsExplicit(w))
            {
                Assert.That(() => AliasProbe.SuffixLinkApart(nav, locus), Throws.TypeOf<ArgumentException>(), at + " (two variables) is no node");
                Assert.That(() => AliasProbe.SuffixLinkAlias(nav, locus), Throws.TypeOf<ArgumentException>(), at + " (one variable) is no node");
                return;
            }
            int[]? link = IndexContractOracle.SuffixLink(w);
            Outcome expected = AliasCatalog.Expect(o, link is not null && o.Children(w).Count > 0 ? link : null);
            Same(AliasProbe.SuffixLinkApart(nav, locus), expected, () => at + " (two variables)");
            Same(AliasProbe.SuffixLinkAlias(nav, locus), expected, () => at + " (one variable)");
        });
    }

    [Test]
    public void GetChildren_TheChildrenMayOverwriteTheVariableThatHoldsTheNode()
    {
        ForEveryLocus((text, o, nav, w, locus) =>
        {
            string at = At(text, w, "GetChildren");
            if (!o.IsExplicit(w))
            {
                Assert.That(() => AliasProbe.ChildrenApart(nav, locus, 3), Throws.TypeOf<ArgumentException>(), at + " (two variables) is no node");
                Assert.That(() => AliasProbe.ChildrenAlias(nav, locus, 3, 1), Throws.TypeOf<ArgumentException>(), at + " (one variable) is no node");
                return;
            }
            IndexLocus[] children = [.. o.Children(w).Select(k => o.LocusOf(k.Child))];
            foreach (int size in new[] { 0, Math.Max(0, children.Length - 1), children.Length, children.Length + 1, children.Length + 3 }.Distinct())
            {
                ChildrenOutcome apart = AliasProbe.ChildrenApart(nav, locus, size);
                if (apart.Count != children.Length) Assert.Fail($"{at} size {size} (two variables): expected {children.Length} children, got {apart.Count}");
                if (!apart.Buffer.SequenceEqual(AliasProbe.ExpectedBuffer(children, size, default, -1))) Assert.Fail($"{at} size {size} (two variables): expected {new ChildrenOutcome(children.Length, AliasProbe.ExpectedBuffer(children, size, default, -1))}, got {apart}");
                for (int slot = 0; slot < size; slot++)
                {
                    ChildrenOutcome alias = AliasProbe.ChildrenAlias(nav, locus, size, slot);
                    IndexLocus[] expected = AliasProbe.ExpectedBuffer(children, size, locus, slot);
                    if (alias.Count != children.Length) Assert.Fail($"{at} size {size} slot {slot} (one variable): expected {children.Length} children, got {alias.Count}");
                    if (!alias.Buffer.SequenceEqual(expected)) Assert.Fail($"{at} size {size} slot {slot} (one variable): expected {new ChildrenOutcome(children.Length, expected)}, got {alias}");
                }
            }
        });
    }

    [Test]
    public void GetChildren_TheChildrenOfTheSecondChild_MayBeWrittenOverTheArrayThatHoldsIt()
    {
        // the form of the review: the children array is reused to descend, kids = children of the root, then the children of kids[1] into kids
        ISuffixTreeNavigation nav = Nav("abracadabra");
        var kids = new IndexLocus[8];
        int count = nav.GetChildren(nav.Root, kids);
        Assert.That(count, Is.EqualTo(OracleOf("abracadabra").Children([]).Count), "the terminator leaf, a, b, c, d, r");
        IndexLocus second = kids[1];
        var expected = new IndexLocus[8];
        int expectedCount = nav.GetChildren(second, expected);
        Assert.That(expectedCount, Is.GreaterThan(2), "the node of 'a' has the children $, b, c, d");
        Assert.That(nav.GetChildren(in kids[1], kids), Is.EqualTo(expectedCount));
        Assert.That(kids.Take(expectedCount), Is.EqualTo(expected.Take(expectedCount)));
    }

    [Test]
    public void Locate_ThePositionsMayBeWrittenOverTheBytesOfTheLocus()
    {
        // not an in/out pair of IndexLocus, but the same hazard: the destination memory and the locus may be one object
        foreach (string text in Texts)
        {
            IndexContractOracle o = OracleOf(text);
            ISuffixTreeNavigation nav = Nav(text);
            foreach (int[] w in AliasCatalog.Loci(o))
            {
                IndexLocus locus = o.LocusOf(w);
                long[] all = [.. o.Positions(w).OrderBy(p => o.Isa[p]).Select(p => (long)p)];
                foreach (int cells in new[] { 1, 2, 4 })
                {
                    foreach (long skip in new long[] { 0, 1 })
                    {
                        IndexLocus[] memory = new IndexLocus[cells];
                        memory[0] = locus;
                        Span<long> positions = MemoryMarshal.Cast<IndexLocus, long>(memory.AsSpan());
                        int written = nav.Locate(in memory[0], skip, positions);
                        long[] expected = [.. all.Skip((int)skip).Take(positions.Length)];
                        Assert.That(written, Is.EqualTo(expected.Length), $"'{text}' {IndexContractOracle.Render(w)} cells {cells} skip {skip}");
                        Assert.That(positions[..written].ToArray(), Is.EqualTo(expected), $"'{text}' {IndexContractOracle.Render(w)} cells {cells} skip {skip}");
                    }
                }
            }
        }
    }
}
